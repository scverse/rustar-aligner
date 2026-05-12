# STAR Suffix Array Implementation Research

**Research Date:** 2026-02-06
**Purpose:** Port STAR's suffix array construction to Rust for rustar-aligner Phase 3

This document provides a detailed analysis of STAR's suffix array construction implementation based on examination of the STAR source code (https://github.com/alexdobin/STAR).

---

## 1. Suffix Array Construction Algorithm

### Algorithm Choice: **Standard qsort with Custom Comparator**

STAR does **NOT** use a specialized suffix array construction algorithm like SA-IS or DC3. Instead, it uses:

- **Standard C library `qsort`** (from `<stdlib.h>`)
- **Custom comparator function** (`funCompareSuffixes`)
- **Prefix-based bucketing** for parallelization and memory efficiency

**Source Location:** `source/Genome_genomeGenerate.cpp`, lines 213-349

### High-Level Algorithm:

1. **Prepare genome**: Reverse the genome array in-place for backward suffix comparison
2. **Calculate 4-nt (2-byte) prefixes** for all suffix positions
3. **Distribute suffixes into chunks** based on prefix buckets to fit in RAM
4. **Sort each chunk** in parallel using `qsort` with `funCompareSuffixes`
5. **Write sorted chunks** to disk as temporary files
6. **Read chunks back** and pack into final PackedArray format

### Key Implementation Details:

```cpp
// From Genome_genomeGenerate.cpp, line 215-217
for (uint ii=0;ii<nGenome;ii++) {
    swap(G[2*nGenome-1-ii],G[ii]);  // Reverse genome for sorting
};
globalG=G;
globalL=pGe.gSuffixLengthMax/sizeof(uint);  // Max comparison length in uint64 words
```

**Why reverse the genome?**
- The comparator works *backwards* from each suffix position (reads `G[ii], G[ii-1], G[ii-2], ...`)
- This allows comparing suffixes as 8-byte chunks using `uint64` comparisons
- More cache-efficient than forward iteration

---

## 2. The Comparator Function: `funCompareSuffixes`

**Source:** `source/Genome_genomeGenerate.cpp`, lines 29-89

### Function Signature:
```cpp
inline int funCompareSuffixes(const void *a, const void *b)
```

### What It Compares:

The comparator receives pointers to **suffix positions** (stored as `uint` indices into the genome array `G`). It compares the genomic sequences starting at those positions.

### Detailed Behavior:

#### 2.1 8-Byte Word Comparisons

```cpp
uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));

uint jj=0;
uint va=0, vb=0;

while (jj < globalL) {
    va=*(ga-jj);  // Read 8-byte word backwards
    vb=*(gb-jj);

    // ... comparison logic ...
    jj++;
};
```

**Key Points:**
- Compares suffixes 8 bytes (64 bits) at a time
- Reads **backwards** from the suffix position (hence `ga-jj`)
- `globalG-7LLU` offset ensures the pointer arithmetic works correctly
- `globalL` = `gSuffixLengthMax / sizeof(uint)` = maximum comparison depth (default 100 bytes = ~12-13 uint64 words)

#### 2.2 Sentinel Detection (Chromosome Boundaries)

```cpp
#define has5(v) ((((v)^0x0505050505050505) - 0x0101010101010101) & ~((v)^0x0505050505050505) & 0x8080808080808080)

if (has5(va) && has5(vb)) {
    // There is a '5' (sentinel/spacer) in the 8-byte word
    // Switch to byte-by-byte comparison
    va1=(uint8*) &va;
    vb1=(uint8*) &vb;
    for (ii=7;ii>=0;ii--) {
        if (va1[ii]>vb1[ii]) return 1;
        else if (va1[ii]<vb1[ii]) return -1;
        else if (va1[ii]==5) {
            // Reached chromosome boundary - break tie using position
            if (*((uint*)a) > *((uint*)b)) return -1;
            else return 1;
        };
    };
}
```

**The `has5` Macro:**
- Brilliant bit-twiddling hack to detect byte value `5` in a 64-bit word
- `5` is used as a **sentinel/spacer character** between chromosomes
- When both suffixes hit `5` at the same byte, they've reached a chromosome boundary
- Tie-breaking: Uses **anti-stable order** (reverse position order) since indices are stored backwards

#### 2.3 Fast Path for No Sentinels

```cpp
else {
    // No sentinel in this 8-byte chunk - simple comparison
    if (va>vb) return 1;
    else if (va<vb) return -1;
};
```

#### 2.4 Final Tie-Breaking

If suffixes are equal up to `globalL` depth:

```cpp
// Anti-stable order since indices are sorted in reverse order
if (*((uint*)a) > *((uint*)b)) return -1;
else return 1;
```

---

## 3. PackedArray Binary Format for SA

**Source Files:**
- `source/PackedArray.h`
- `source/PackedArray.cpp`

### Purpose
Store the suffix array using **variable bits per entry** to save memory (e.g., 33 bits instead of 64 bits for human genome).

### Key Fields:

```cpp
class PackedArray {
    uint wordLength;        // Bits per element
    uint length;            // Number of elements
    uint lengthByte;        // Total bytes allocated
    char* charArray;        // Actual data storage
    uint bitRecMask;        // Mask for extracting values
    uint wordCompLength;    // Complement of wordLength (for shifting)
};
```

### Bits Per Entry Formula

```cpp
void PackedArray::defineBits(uint Nbits, uint lengthIn) {
    wordLength = Nbits;                              // Store bits per element
    wordCompLength = sizeof(uint)*8LLU - wordLength; // = 64 - wordLength
    bitRecMask = (~0LLU) >> wordCompLength;          // Mask to extract value
    length = lengthIn;                                // Number of elements
    lengthByte = (length-1)*wordLength/8LLU + sizeof(uint);
}
```

**For the suffix array:**
```cpp
// From Genome_genomeGenerate.cpp, line 196
SA.defineBits(GstrandBit+1, nSA);
```

Where:
```cpp
// Line 191
GstrandBit = (char)(uint)floor(log(nGenome + P.limitSjdbInsertNsj*sjdbLength)/log(2)) + 1;
if (GstrandBit<32) GstrandBit=32;
```

**Example for human genome:**
- `nGenome ≈ 3.1e9 bases`
- `GstrandBit = floor(log2(3.1e9)) + 1 = 32`
- `wordLength = GstrandBit + 1 = 33 bits per SA entry`

### Write Operation

```cpp
void PackedArray::writePacked(uint jj, uint x) {
    uint b = jj * wordLength;        // Bit position
    uint B = b / 8LLU;                // Byte offset
    uint S = b % 8LLU;                // Bit shift within byte

    x = x << S;                       // Shift value to position
    uint* a1 = (uint*) (charArray+B); // Get 64-bit word at position
    *a1 = ((*a1) & ~(bitRecMask<<S)) | x;  // Clear bits and write value
}
```

### Read Operation

```cpp
inline uint PackedArray::operator[](uint ii) {
    uint b = ii * wordLength;         // Bit position
    uint B = b / 8;                   // Byte offset
    uint S = b % 8;                   // Bit shift

    uint a1 = *((uint*)(charArray+B));
    a1 = ((a1>>S)<<wordCompLength)>>wordCompLength;  // Extract and mask
    return a1;
}
```

### Byte Layout and Endianness

- **Endianness:** Uses **native endianness** (little-endian on x86-64)
- **Bit packing:** Values can span byte boundaries
- **Word alignment:** Allocates extra `sizeof(uint)` at the end to avoid overruns
- **Storage format:** Raw binary dump of `charArray` to disk

---

## 4. Prefix Bucketing for Parallelism

**Source:** `source/Genome_genomeGenerate.cpp`, lines 220-256

### Prefix Length: **4 nucleotides = 2 bytes = 16 bits**

```cpp
uint indPrefN = 1LLU << 16;  // 65,536 possible 4-nt prefixes
uint* indPrefCount = new uint[indPrefN];
```

### Bucketing Strategy:

```cpp
for (uint ii=0; ii<2*nGenome; ii+=pGe.gSAsparseD) {
    if (G[ii]<4) {  // Valid nucleotide (not N)
        // Extract 4-nt prefix (reading backwards)
        uint p1 = (G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
        indPrefCount[p1]++;
        nSA++;
    };
};
```

**Prefix encoding:** Each nucleotide is 2 bits: `ACGT = 0,1,2,3`

### Dynamic Chunk Sizing:

```cpp
uint saChunkSize = (P.limitGenomeGenerateRAM - nG1alloc) / 8 / P.runThreadN;
saChunkSize = saChunkSize * 6 / 10;  // Allow 40% overhead for qsort
```

- Calculates how many SA indices fit in available RAM
- Accounts for `qsort` temporary memory (~40% overhead)
- Divides work across `P.runThreadN` threads

### Chunk Assignment:

```cpp
indPrefStart[0] = 0;
saChunkN = 0;
uint chunkSize1 = indPrefCount[0];

for (uint ii=1; ii<indPrefN; ii++) {
    chunkSize1 += indPrefCount[ii];
    if (chunkSize1 > saChunkSize) {
        saChunkN++;
        indPrefStart[saChunkN] = ii;
        indPrefChunkCount[saChunkN-1] = chunkSize1 - indPrefCount[ii];
        chunkSize1 = indPrefCount[ii];
    };
};
```

**Result:** Each chunk contains suffixes with consecutive prefix values, sized to fit in RAM.

### Parallel Sorting:

```cpp
#pragma omp parallel for num_threads(P.runThreadN) ordered schedule(dynamic,1)
for (int iChunk=0; iChunk < (int)saChunkN; iChunk++) {
    uint* saChunk = new uint[indPrefChunkCount[iChunk]];

    // Fill chunk with matching suffixes
    for (uint ii=0,jj=0; ii<2*nGenome; ii+=pGe.gSAsparseD) {
        if (G[ii]<4) {
            uint p1 = (G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
            if (p1>=indPrefStart[iChunk] && p1<indPrefStart[iChunk+1]) {
                saChunk[jj] = ii;
                jj++;
            };
        };
    };

    // Sort the chunk
    qsort(saChunk, indPrefChunkCount[iChunk], sizeof(saChunk[0]), funCompareSuffixes);

    // Convert back to forward indexing
    for (uint ii=0; ii<indPrefChunkCount[iChunk]; ii++) {
        saChunk[ii] = 2*nGenome - 1 - saChunk[ii];
    };

    // Write to disk
    string chunkFileName = pGe.gDir + "/SA_" + to_string((uint)iChunk);
    ofstream saChunkFile(chunkFileName);
    fstreamWriteBig(saChunkFile, (char*)saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk], ...);
    saChunkFile.close();
    delete[] saChunk;
};
```

---

## 5. SA Index (SAindex) Structure

**Source Files:**
- `source/genomeSAindex.cpp`
- `source/genomeSAindex.h`

### Purpose
Create a **lookup table** that maps k-mer prefixes to their starting positions in the sorted suffix array, enabling fast binary search during alignment.

### SAindex Parameters

```cpp
uint gSAindexNbases;  // Default: 14 (14-mer prefixes)
```

From the STAR manual:
- Typical range: 10-15
- Longer = more memory, faster search
- Shorter = less memory, slower search

### SAindex Size Calculation

```cpp
mapGen.genomeSAindexStart = new uint[mapGen.pGe.gSAindexNbases+1];
mapGen.genomeSAindexStart[0] = 0;

for (uint ii=1; ii<=mapGen.pGe.gSAindexNbases; ii++) {
    mapGen.genomeSAindexStart[ii] = mapGen.genomeSAindexStart[ii-1] + (1LLU<<(2*ii));
};

mapGen.nSAi = mapGen.genomeSAindexStart[mapGen.pGe.gSAindexNbases];
```

**For `gSAindexNbases=14`:**
- 1-mers: 4^1 = 4 entries
- 2-mers: 4^2 = 16 entries
- 3-mers: 4^3 = 64 entries
- ...
- 14-mers: 4^14 = 268,435,456 entries
- **Total:** `nSAi = (4^15 - 4) / 3 ≈ 357 million entries`

### SAindex Entry Format

Each SAindex entry stores the **first SA position** where a k-mer prefix appears.

**Bits per SAindex entry:**
```cpp
SAi.defineBits(mapGen.GstrandBit+3, mapGen.nSAi);
```

- `GstrandBit + 1`: SA position
- `+1 bit`: "Contains N" flag (`SAiMarkNbit`)
- `+1 bit`: "Prefix absent" flag (`SAiMarkAbsentBit`)
- **Total: 35 bits per entry** (for human genome)

### Bit Flags:

```cpp
// Line 28-34 in genomeSAindex.cpp
mapGen.SAiMarkNbit = mapGen.GstrandBit + 1;
mapGen.SAiMarkAbsentBit = mapGen.GstrandBit + 2;

mapGen.SAiMarkNmaskC = 1LLU << mapGen.SAiMarkNbit;           // Set flag
mapGen.SAiMarkNmask = ~mapGen.SAiMarkNmaskC;                 // Clear flag
mapGen.SAiMarkAbsentMaskC = 1LLU << mapGen.SAiMarkAbsentBit; // Set flag
mapGen.SAiMarkAbsentMask = ~mapGen.SAiMarkAbsentMaskC;       // Clear flag
```

### SAindex Generation Algorithm

**Main function:** `genomeSAindexChunk()` in `genomeSAindex.cpp`, lines 117-176

```cpp
void genomeSAindexChunk(char * G, PackedArray & SA, Parameters & P,
                        PackedArray & SAi, uint iSA1, uint iSA2, Genome &mapGen)
{
    uint* ind0 = new uint[mapGen.pGe.gSAindexNbases];

    // Initialize to -1 (in case some prefixes don't exist)
    for (uint ii=0; ii<mapGen.pGe.gSAindexNbases; ii++) {
        ind0[ii] = -1;
    };

    uint isaStep = mapGen.nSA / (1llu<<(2*mapGen.pGe.gSAindexNbases)) + 1;
    uint isa = iSA1;
    int iL4;

    // Calculate initial k-mer prefix
    uint indFull = funCalcSAiFromSA(G, SA, mapGen, isa, mapGen.pGe.gSAindexNbases, iL4);

    while (isa <= iSA2) {
        // For each prefix length (1-mer, 2-mer, ..., k-mer)
        for (uint iL=0; iL<mapGen.pGe.gSAindexNbases; iL++) {
            // Extract (iL+1)-mer prefix
            uint indPref = indFull >> (2*(mapGen.pGe.gSAindexNbases-1-iL));

            if ((int)iL == iL4) {
                // This suffix contains N - mark all longer prefixes
                for (uint iL1=iL; iL1<mapGen.pGe.gSAindexNbases; iL1++) {
                    SAi.writePacked(
                        mapGen.genomeSAindexStart[iL1] + ind0[iL1],
                        SAi[mapGen.genomeSAindexStart[iL1] + ind0[iL1]] | mapGen.SAiMarkNmaskC
                    );
                };
                break;
            };

            if (indPref > ind0[iL] || isa==0) {
                // New prefix - record SA position
                SAi.writePacked(mapGen.genomeSAindexStart[iL] + indPref, isa);

                // Fill gaps (absent prefixes) with current position marked as absent
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {
                    SAi.writePacked(
                        mapGen.genomeSAindexStart[iL] + ii,
                        isa | mapGen.SAiMarkAbsentMaskC
                    );
                };
                ind0[iL] = indPref;
            }
        };

        // Find next different suffix
        funSAiFindNextIndex(G, SA, isaStep, isa, indFull, iL4, mapGen);
    };

    // Fill remaining absent prefixes at the end
    for (uint iL=0; iL<mapGen.pGe.gSAindexNbases; iL++) {
        for (uint ii=mapGen.genomeSAindexStart[iL]+ind0[iL]+1;
             ii<mapGen.genomeSAindexStart[iL+1]; ii++) {
            SAi.writePacked(ii, mapGen.nSA | mapGen.SAiMarkAbsentMaskC);
        };
    };

    delete[] ind0;
}
```

### Helper: `funCalcSAiFromSA()`

**Source:** `SuffixArrayFuns.cpp`, lines 353-395

Extracts a k-mer from the genome at a given SA position:

```cpp
uint funCalcSAiFromSA(char* gSeq, PackedArray& gSA, Genome &mapGen,
                      uint iSA, int L, int & iL4)
{
    uint SAstr = gSA[iSA];
    bool dirG = (SAstr>>mapGen.GstrandBit) == 0;  // Forward/reverse strand
    SAstr &= mapGen.GstrandMask;                   // Mask off strand bit

    iL4 = -1;
    register uint saind = 0;

    if (dirG) {
        // Forward strand
        register uint128 g1 = *((uint128*)(gSeq+SAstr));
        for (int ii=0; ii<L; ii++) {
            register char g2 = (char)g1;
            if (g2 > 3) {  // Hit N or sentinel
                iL4 = ii;
                saind <<= 2*(L-ii);
                return saind;
            };
            saind = saind<<2;
            saind += g2;
            g1 = g1>>8;
        };
        return saind;
    } else {
        // Reverse strand - complement and reverse
        register uint128 g1 = *((uint128*)(gSeq+mapGen.nGenome-SAstr-16));
        for (int ii=0; ii<L; ii++) {
            register char g2 = (char)(g1>>(8*(15-ii)));
            if (g2 > 3) {
                iL4 = ii;
                saind <<= 2*(L-ii);
                return saind;
            };
            saind = saind<<2;
            saind += 3-g2;  // Complement
        };
        return saind;
    };
}
```

**Key points:**
- Reads 16 bytes at once using `uint128`
- Returns k-mer as packed integer (2 bits per base)
- Sets `iL4` to position of first N/sentinel (or -1 if none)

### SAindex Binary File Format

**File:** `{genomeDir}/SAindex`

**Write code:** `Genome_genomeGenerate.cpp`, lines 398-404

```cpp
ofstream & SAiOut = ofstrOpen(pGe.gDir+"/SAindex", ERROR_OUT, P);

// Header: number of bases in index
fstreamWriteBig(SAiOut, (char*)&pGe.gSAindexNbases,
                sizeof(pGe.gSAindexNbases), ...);

// Prefix start positions array (length gSAindexNbases+1)
fstreamWriteBig(SAiOut, (char*)genomeSAindexStart,
                sizeof(genomeSAindexStart[0])*(pGe.gSAindexNbases+1), ...);

// Packed SAindex data
fstreamWriteBig(SAiOut, SAi.charArray, SAi.lengthByte, ...);

SAiOut.close();
```

**Binary layout:**
```
[uint64: gSAindexNbases]                     // 8 bytes
[uint64[gSAindexNbases+1]: genomeSAindexStart]  // 8*(gSAindexNbases+1) bytes
[byte[SAi.lengthByte]: packed SAi data]      // Variable length
```

---

## 6. GstrandBit Parameter

### Purpose
Encodes both **genome position** and **strand** in a single integer.

### Bit Layout:

For a 33-bit SA entry (human genome):
```
Bit 32: Strand bit (0=forward, 1=reverse)
Bits 0-31: Genome position
```

### Calculation:

```cpp
// Line 191 in Genome_genomeGenerate.cpp
GstrandBit = (char)(uint)floor(log(nGenome + P.limitSjdbInsertNsj*sjdbLength)/log(2)) + 1;
if (GstrandBit < 32) GstrandBit = 32;

// Line 195
GstrandMask = ~(1LLU << GstrandBit);
```

**Example:**
- Human genome: ~3.1 billion bases
- `log2(3.1e9) ≈ 31.5`
- `GstrandBit = 32`
- `GstrandMask = 0xFFFFFFFF` (mask off bit 32)

### Encoding/Decoding:

```cpp
// Encode (line 316 in Genome_genomeGenerate.cpp)
SA.writePacked(packedInd+ii,
    (saIn[ii]<nGenome) ? saIn[ii] : ((saIn[ii]-nGenome) | N2bit));

where:
    uint N2bit = 1LLU << GstrandBit;

// Decode (from SuffixArrayFuns.cpp, line 20-22)
uint SAstr = mapGen.SA[iSA];
bool dirG = (SAstr>>mapGen.GstrandBit) == 0;  // Forward = 0, Reverse = 1
SAstr &= mapGen.GstrandMask;                   // Extract position
```

### Genome Organization:

STAR stores the genome as:
```
[0 ... nGenome-1]: Forward strand
[nGenome ... 2*nGenome-1]: Reverse complement
```

During SA construction:
- Suffixes from `[0, nGenome)` are encoded with strand bit = 0
- Suffixes from `[nGenome, 2*nGenome)` are encoded with strand bit = 1 and position adjusted

---

## 7. Practical Details

### 7.1 Memory Usage During SA Construction

**Peak memory:**
```cpp
// Genome storage (forward + reverse)
uint64 genomeBytes = 2 * nGenome;

// Chunk size per thread
uint saChunkSize = (P.limitGenomeGenerateRAM - nG1alloc) / 8 / P.runThreadN * 6 / 10;

// Total memory:
// - Genome: 2*nGenome bytes
// - Per-thread chunk: saChunkSize * 4 bytes (uint)
// - qsort overhead: ~40% extra
```

**For human genome with 32GB RAM:**
- Genome: ~6 GB (2 * 3.1 GB)
- Available for SA: ~26 GB
- With 8 threads: ~3.25 GB per thread
- Chunk size: ~1.95 GB = ~488M SA indices per chunk

### 7.2 Sparsity Parameter

```cpp
uint gSAsparseD;  // Default: 1 (no sparsity)
```

STAR can construct a **sparse suffix array** by only indexing every `gSAsparseD`-th position:

```cpp
for (uint ii=0; ii<2*nGenome; ii+=pGe.gSAsparseD) {
    if (G[ii]<4) {
        nSA++;
    };
};
```

**Trade-offs:**
- `gSAsparseD=1`: Full SA, maximum sensitivity, large memory
- `gSAsparseD=2`: Half the memory, slightly reduced sensitivity
- Typical: 1-2

### 7.3 Sentinel Values

**Sentinel character: `5` (defined as `GENOME_spacingChar`)**

Inserted between chromosomes during genome concatenation:

```cpp
// From genomeScanFastaFiles (not shown in excerpts)
// Chromosomes are concatenated with spacers of 100x char(5)
```

**Why 5?**
- Nucleotides: `A=0, C=1, G=2, T=3, N=4`
- `5` is > all valid bases
- Easy to detect in comparisons
- Used to stop suffix comparisons at chromosome boundaries

### 7.4 Handling Forward/Reverse Strands

**During construction:**

1. **Genome preparation** (line 178-182):
```cpp
for (uint ii=0; ii<nGenome; ii++) {
    G[2*nGenome-1-ii] = G[ii]<4 ? 3-G[ii] : G[ii];  // Reverse complement
};
```

2. **SA indexing** spans both strands:
```cpp
for (uint ii=0; ii<2*nGenome; ii+=pGe.gSAsparseD) {
    if (G[ii]<4) {
        // Add to SA
    };
};
```

3. **Encoding in SA** (line 316):
```cpp
SA.writePacked(packedInd+ii,
    (saIn[ii]<nGenome) ? saIn[ii] : ((saIn[ii]-nGenome) | N2bit));
```

**Result:** SA contains suffixes from both strands, sorted together, with strand bit encoded.

---

## 8. Implementation Roadmap for Rust

### Phase 1: Data Structures

1. **PackedArray implementation:**
   - Variable bit-width storage
   - Read/write methods with bit manipulation
   - Allocation/deallocation

2. **Genome structure:**
   - Forward + reverse complement storage
   - Chromosome boundaries with sentinels
   - GstrandBit calculation

### Phase 2: Comparator

1. **Suffix comparator function:**
   - 8-byte word comparison (use `u64`)
   - Sentinel detection with bit tricks (`has5` macro)
   - Byte-by-byte fallback
   - Anti-stable tie-breaking

2. **Alternative:** Port `funCompareUintAndSuffixesMemcmp` for simpler implementation

### Phase 3: SA Construction

1. **Prefix bucketing:**
   - 4-nt prefix extraction
   - Chunk size calculation based on available RAM
   - Distribute suffixes to chunks

2. **Parallel sorting:**
   - Use `rayon` for parallelism
   - Sort each chunk with custom comparator
   - Write to temporary files

3. **Merging:**
   - Read chunks sequentially
   - Pack into PackedArray format
   - Write final SA to disk

### Phase 4: SAindex Generation

1. **Calculate k-mer prefixes** from SA positions
2. **Build SAindex lookup table:**
   - Track first occurrence of each k-mer
   - Mark absent prefixes
   - Mark prefixes containing N

3. **Write SAindex file:**
   - Header with `gSAindexNbases`
   - Prefix start positions
   - Packed SAi data

### Phase 5: I/O

1. **Binary file writing:**
   - Genome file
   - SA file (packed format)
   - SAindex file (packed format)

2. **Binary file reading:**
   - Load genome
   - Load SA
   - Load SAindex

---

## 9. Key Formulas Reference

### GstrandBit:
```
GstrandBit = max(32, floor(log2(nGenome + limitSjdbInsertNsj*sjdbLength)) + 1)
```

### SA wordLength:
```
SA.wordLength = GstrandBit + 1
```

### SAindex size:
```
nSAi = sum(4^i for i=1 to gSAindexNbases)
     = (4^(gSAindexNbases+1) - 4) / 3
```

### SAindex wordLength:
```
SAi.wordLength = GstrandBit + 3  // +1 for SA index, +2 for flags
```

### PackedArray lengthByte:
```
lengthByte = (length-1) * wordLength / 8 + sizeof(uint)
```

### Chunk size calculation:
```
saChunkSize = (limitGenomeGenerateRAM - genomeBytes) / 8 / numThreads * 0.6
```

---

## 10. Sources

All information derived from the STAR source code repository:
- Repository: https://github.com/alexdobin/STAR
- Version examined: master branch (as of 2026-02-06)

Key files analyzed:
- `source/Genome_genomeGenerate.cpp` - Main SA generation logic
- `source/SuffixArrayFuns.cpp` - SA utility functions
- `source/PackedArray.h` / `.cpp` - Bit-packed array storage
- `source/genomeSAindex.cpp` / `.h` - SAindex generation
- `source/funCompareUintAndSuffixes.cpp` - SA comparator (alternative)
- `source/Genome.h` - Core data structures

---

## Appendix: Complete Comparator Pseudocode

```rust
fn compare_suffixes(a_pos: u32, b_pos: u32, genome: &[u8], max_depth: usize) -> Ordering {
    // Convert to u64 pointers for 8-byte comparisons
    let mut depth = 0;

    while depth < max_depth {
        let offset = depth * 8;

        // Read 8-byte words backwards from suffix positions
        let word_a = read_u64_backwards(&genome, a_pos, offset);
        let word_b = read_u64_backwards(&genome, b_pos, offset);

        // Check if either word contains sentinel (value 5)
        if has_sentinel(word_a) && has_sentinel(word_b) {
            // Byte-by-byte comparison needed
            for byte_idx in (0..8).rev() {
                let byte_a = ((word_a >> (byte_idx * 8)) & 0xFF) as u8;
                let byte_b = ((word_b >> (byte_idx * 8)) & 0xFF) as u8;

                match byte_a.cmp(&byte_b) {
                    Ordering::Equal if byte_a == 5 => {
                        // Reached chromosome boundary - anti-stable sort
                        return b_pos.cmp(&a_pos);  // Note: reversed
                    }
                    Ordering::Equal => continue,
                    other => return other,
                }
            }
        } else {
            // Fast path: no sentinel, compare as u64
            match word_a.cmp(&word_b) {
                Ordering::Equal => {
                    depth += 1;
                    continue;
                }
                other => return other,
            }
        }

        depth += 1;
    }

    // Suffixes equal up to max_depth - tie-break with position (anti-stable)
    b_pos.cmp(&a_pos)
}

fn has_sentinel(word: u64) -> bool {
    // Detect byte value 5 in u64 word using bit tricks
    const SENTINEL_PATTERN: u64 = 0x0505050505050505;
    const SUB_MASK: u64 = 0x0101010101010101;
    const DETECT_MASK: u64 = 0x8080808080808080;

    let xor = word ^ SENTINEL_PATTERN;
    let sub = xor.wrapping_sub(SUB_MASK);
    let detect = (!xor) & sub & DETECT_MASK;

    detect != 0
}
```

---

**End of Research Document**
