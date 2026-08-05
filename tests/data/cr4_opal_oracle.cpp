// CR4 TSO-clip oracle: STAR's real Opal call, transcribed from ClipCR4.cpp +
// ClipMate_clipChunk.cpp, emitting `read \t clip \t score \t L` per read.
//
// Generates a deterministic, boundary-heavy CR4 read corpus, runs it through
// Opal in STAR's 64-read chunks with STAR's exact parameters, and applies
// STAR's clip gate. The TSV it prints is the committed expected output for
// rustar-aligner's `cr4_tso_matches_star_opal_oracle` test.
//
// Regenerate (needs a STAR checkout for opal.cpp — MIT, (c) Alexander Dobin):
//
//   OPAL=/path/to/STAR/source/opal
//   g++ -O2 -std=c++11 -I"$OPAL" tests/data/cr4_opal_oracle.cpp "$OPAL/opal.cpp" \
//       -o /tmp/cr4_opal_oracle
//   /tmp/cr4_opal_oracle > tests/data/cr4_opal_oracle.tsv
//
// The corpus is driven by a fixed-seed LCG, so the output is reproducible; any
// change to the generator invalidates the committed TSV and must be regenerated.
#include "opal.h"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

static const char *TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG";
static const uint32_t READ_LEN = 91;
static const int DBN = 64;

static uint64_t st = 0x51ED2701ULL;
static uint64_t lcg() {
    st = st * 6364136223846793005ULL + 1442695040888963407ULL;
    return st >> 33;
}
static char base(uint64_t v) { return "ACGT"[v % 4]; }

// A corpus that straddles STAR's S<20 / S==20&&L>26 / S==21&&L>30 decision
// boundary: TSO with 0..14 mismatches, 5' shifts, truncations (overlap mode's
// whole point), indels, embedded Ns, polyA, and pure random.
static std::vector<std::string> corpus() {
    std::vector<std::string> reads;
    std::string tso(TSO);

    for (int i = 0; i < 600; i++) {
        std::string r = tso;
        int k = lcg() % 15;                       // 0..14 mismatches
        for (int j = 0; j < k; j++) {
            size_t p = lcg() % r.size();
            r[p] = base(lcg());
        }
        if (lcg() % 4 == 0) {                     // truncate the TSO (partial overlap)
            r = r.substr(lcg() % 12);
        }
        if (lcg() % 5 == 0) {                     // indel inside the adapter
            size_t p = lcg() % r.size();
            if (lcg() % 2) r.erase(p, 1); else r.insert(p, 1, base(lcg()));
        }
        size_t shift = lcg() % 6;                 // 5' offset before the TSO
        std::string read;
        for (size_t j = 0; j < shift; j++) read += base(lcg());
        read += r;
        while (read.size() < READ_LEN) read += base(lcg());
        read = read.substr(0, READ_LEN);
        if (lcg() % 8 == 0) read[lcg() % read.size()] = 'N';
        reads.push_back(read);
    }
    for (int i = 0; i < 150; i++) {               // pure random
        std::string read;
        for (uint32_t j = 0; j < READ_LEN; j++) read += base(lcg());
        reads.push_back(read);
    }
    for (int i = 0; i < 50; i++) {                // polyA / polyT / N-heavy
        char c = (i % 3 == 0) ? 'A' : (i % 3 == 1) ? 'T' : 'N';
        reads.push_back(std::string(READ_LEN, c));
    }
    for (int i = 0; i < 100; i++) {               // short reads (STAR N-pads to 91)
        size_t len = 5 + (lcg() % 80);
        std::string read = tso.substr(0, std::min(len, tso.size()));
        while (read.size() < len) read += base(lcg());
        reads.push_back(read);
    }

    // ---- poly-A corpus, for polyTail3p (appended last so the reads above keep
    // their sequences and the TSO expectations stay stable) ----
    // Straddles polyTail3p's own boundaries: the score>=20 floor, the 70%
    // density rule, the ib-score>27 give-up, and the seqLen<20 early return.
    for (size_t tail : {0, 1, 5, 9, 10, 15, 19, 20, 21, 25, 30, 40, 60, 91}) {
        std::string head;
        while (head.size() + tail < READ_LEN) head += "CGT"[head.size() % 3];
        reads.push_back(head + std::string(tail, 'A'));
    }
    // A-rich upstream: the scan does not stop at the tail boundary.
    for (size_t tail : {20, 25, 30}) {
        std::string head;
        while (head.size() + tail < READ_LEN) head += "ACGT"[head.size() % 4];
        reads.push_back(head + std::string(tail, 'A'));
    }
    // Interrupted tails: one or more non-A inside an otherwise clean run.
    for (size_t gap = 1; gap <= 12; gap++) {
        std::string read(READ_LEN - 40, 'C');
        std::string tail(40, 'A');
        for (size_t j = 0; j < gap && j * 3 + 1 < tail.size(); j++) tail[j * 3 + 1] = 'C';
        reads.push_back(read + tail);
    }
    // Short reads around the seqLen<20 floor, and TSO+polyA together.
    for (size_t len : {5, 15, 19, 20, 21, 25}) reads.push_back(std::string(len, 'A'));
    for (size_t tail : {10, 25, 40}) {
        std::string read = tso;
        while (read.size() + tail < READ_LEN) read += base(lcg());
        reads.push_back(read + std::string(tail, 'A'));
    }
    return reads;
}

/// STAR `ClipCR4::polyTail3p`, transcribed verbatim. `seq` is numeric (A == 0)
/// and is the read *after* the 5' TSO clip, matching `ClipMate::clip`'s order
/// (type 10 shifts the sequence and shrinks Lread, then type 11 runs on it).
static uint32_t polyTail3p(const uint8_t *seq, uint32_t seqLen) {
    if (seqLen < 20) return 0;
    uint32_t ib1 = seqLen - 1;
    int32_t score = 0, score1 = 0;
    for (uint32_t ib = 1; ib <= seqLen; ib++) {
        if (seq[seqLen - ib] == 0) {
            score += 1;
            if (score * 10 >= (int)ib * 7) { ib1 = ib; score1 = score; }
        } else {
            score -= 2;
            if (ib - score > 27) break;
        }
    }
    if (score1 < 20) ib1 = 0;
    return ib1;
}

static void encode(const std::string &s, uint8_t *dst) {
    uint32_t minLen = std::min((uint32_t)s.size(), READ_LEN);
    for (uint32_t i = 0; i < minLen; i++) {
        switch (s[i]) {
            case 'A': dst[i] = 0; break;
            case 'C': dst[i] = 1; break;
            case 'G': dst[i] = 2; break;
            case 'T': dst[i] = 3; break;
            default:  dst[i] = 4;
        }
    }
    if (s.size() < READ_LEN) memset(dst + s.size(), 4, READ_LEN - s.size());
}

int main() {
    // STAR's ClipCR4 constructor, verbatim.
    std::vector<int> scoreMatrix = {
         1, -2, -2, -2, -2,
        -2,  1, -2, -2, -2,
        -2, -2,  1, -2, -2,
        -2, -2, -2,  1, -2,
        -2, -2, -2, -2,  0
    };
    int alphabetLength = 5, gapOpen = 2, gapExt = 2;

    std::vector<uint8_t> dbSeqArr(DBN * READ_LEN);
    std::vector<uint8_t *> dbSeqs(DBN);
    std::vector<int> dbSeqsLen(DBN, (int)READ_LEN);
    for (int i = 0; i < DBN; i++) dbSeqs[i] = dbSeqArr.data() + i * READ_LEN;

    std::vector<OpalSearchResult> opalRes(DBN);
    std::vector<OpalSearchResult *> opalResP(DBN);
    for (int i = 0; i < DBN; i++) opalResP[i] = &opalRes[i];

    std::vector<uint8_t> query;
    for (const char *p = TSO; *p; p++) {
        query.push_back(*p == 'A' ? 0 : *p == 'C' ? 1 : *p == 'G' ? 2 : *p == 'T' ? 3 : 4);
    }

    std::vector<std::string> reads = corpus();

    for (size_t off = 0; off < reads.size(); off += DBN) {
        int dbN1 = (int)std::min((size_t)DBN, reads.size() - off);
        for (int i = 0; i < dbN1; i++) encode(reads[off + i], dbSeqs[i]);
        for (int i = 0; i < dbN1; i++) opalInitSearchResult(opalResP[i]);

        opalSearchDatabase(query.data(), (int)query.size(), dbSeqs.data(), dbN1,
                           dbSeqsLen.data(), gapOpen, gapExt, scoreMatrix.data(),
                           alphabetLength, opalResP.data(),
                           OPAL_SEARCH_SCORE_END, OPAL_MODE_OV, OPAL_OVERFLOW_BUCKETS);

        for (int i = 0; i < dbN1; i++) {
            int L = opalRes[i].endLocationTarget + 1;
            int S = opalRes[i].score;
            bool L0 = S < 20 || (S == 20 && L > 26) || (S == 21 && L > 30);
            // STAR then applies min(clippedInfo, Lread) in ClipMate_clip.cpp.
            const std::string &r = reads[off + i];
            uint32_t clip = L0 ? 0 : (uint32_t)L;
            clip = std::min(clip, (uint32_t)r.size());

            // 3' poly-A, on the read as it stands after the 5' clip.
            std::vector<uint8_t> enc(r.size());
            for (size_t j = 0; j < r.size(); j++) {
                enc[j] = r[j] == 'A' ? 0 : r[j] == 'C' ? 1 : r[j] == 'G' ? 2 : r[j] == 'T' ? 3 : 4;
            }
            uint32_t polya = polyTail3p(enc.data() + clip, (uint32_t)(r.size() - clip));

            printf("%s\t%u\t%d\t%d\t%u\n", r.c_str(), clip, S, L, polya);
        }
    }
    return 0;
}
