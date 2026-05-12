# Testing Framework Quick Start

Get started with the rustar-aligner testing framework in 3 steps.

## Step 1: Verify Setup

```bash
cd test
./verify_framework.sh
```

**Expected output:**
```
✓ All checks passed (17/17)
✓ Test framework is ready!
```

If any checks fail, follow the suggestions to fix prerequisites.

---

## Step 2: Run Your First Test

```bash
./run_tests.sh yeast_100
```

**What happens:**
1. Runs STAR on 100 yeast reads → `results/TIMESTAMP_yeast_100/star/`
2. Runs rustar-aligner on same reads → `results/TIMESTAMP_yeast_100/rustar/`
3. Compares outputs → `results/TIMESTAMP_yeast_100/comparison/`
4. Prints summary report

**Example output:**
```
[23:45:12] Running test: yeast_100
[23:45:15] Running STAR...
[23:45:18] STAR completed successfully
[23:45:18] Running rustar-aligner...
[23:45:20] rustar-aligner completed successfully
[23:45:20] Comparing outputs...

=== Statistics Comparison: rustar-aligner vs STAR ===
  Total reads                      100 vs      100  ✓
  Mapped reads                      91 vs       92  ⚠ (+1)
  Unique mapped                     82 vs       82  ✓
  Multi-mapped                       9 vs       10  ⚠ (+1)
  Mapping rate                    91.00% vs  92.00%

=== Record-Level Comparison ===
  Common reads: 100
  Summary: 99/100 reads match (99.0%)

Status: PASS ✓
Alignment outputs match within 1.0% tolerance
```

---

## Step 3: Save Golden Outputs

If the test passed and you're satisfied with the results:

```bash
./save_golden.sh yeast_100
```

This saves the rustar-aligner output as a golden reference for future regression testing.

**Files created:**
- `golden/yeast_100/stats.json` — Alignment statistics
- `golden/yeast_100/junctions.json` — Junction data
- `golden/yeast_100/metadata.json` — Test metadata

**Version control:**
```bash
git add golden/yeast_100/
git commit -m "Add yeast_100 golden output"
```

---

## Common Workflows

### Run All Tests

```bash
./run_tests.sh --all
```

Runs all 6 test cases (100, 1K, 1K paired, 1K two-pass, 1K BAM, 10K).

### Run Multiple Tests

```bash
./run_tests.sh yeast_100 yeast_1k yeast_10k
```

### Run Tests in Parallel

```bash
./run_tests.sh --all --parallel
```

Faster! Runs independent tests concurrently.

### Run CI Tests (Fast)

```bash
./ci.sh
```

Quick validation for pre-commit or CI/CD:
1. Build rustar-aligner
2. Run unit tests
3. Run yeast_100 integration test
4. Compare against golden outputs

Exit code: 0=pass, 1=regression, 2=build fail, 3=execution fail

---

## When Tests Fail

### 1. View the Comparison Report

```bash
cat results/TIMESTAMP_yeast_1k/comparison/summary.txt
```

### 2. Investigate Discrepancies

```bash
./investigate.sh results/TIMESTAMP_yeast_1k/comparison/alignment_diff.txt
```

**Output:**
- Discrepant read names and issues
- Likely causes and code locations
- Suggestions for debugging

### 3. Debug with Logging

```bash
RUST_LOG=debug ../target/release/rustar-aligner \
  --runMode alignReads \
  --genomeDir data/small/yeast/indices \
  --readFilesIn debug/discrepant_reads.fq \
  --outFileNamePrefix debug/ \
  --outSAMtype SAM
```

### 4. Compare STAR Source

Check `STAR_REFERENCE.md` for:
- rustar-aligner → STAR module mapping
- Algorithm descriptions
- Known differences

---

## Test Cases Reference

| Name | Reads | Time | Purpose |
|------|-------|------|---------|
| `yeast_100` | 100 | ~5s | Fast smoke test |
| `yeast_1k` | 1,000 | ~10s | Standard validation |
| `yeast_1k_bam` | 1,000 | ~10s | BAM output test |
| `yeast_1k_paired` | 1,000 | ~10s | Paired-end test |
| `yeast_1k_twopass` | 1,000 | ~15s | Two-pass mode test |
| `yeast_10k` | 10,000 | ~30s | Performance test |

---

## Tips

### Keep Only Recent Results

```bash
# Results can be large - clean up old runs
rm -rf results/2026020*  # Delete old timestamps
```

### Compare Specific Files Manually

```bash
python compare_sam.py \
  --star star_output/Aligned.out.sam \
  --rustar-aligner rustar_aligner_output/Aligned.out.sam \
  --tolerance 0.01 \
  --verbose
```

### Check Prerequisites Anytime

```bash
./verify_framework.sh
```

### Get Help

```bash
./run_tests.sh --help
```

---

## Expected Results

**Typical results for yeast_1k test:**
- Total reads: 1000
- Mapped: ~900 (90%)
- Unique: ~820 (82% of mapped)
- Multi-mapped: ~80 (8% of mapped)
- Unmapped: ~100 (10%)

**Pass criteria:**
- ≥99% of reads have matching alignments
- Statistics within 1% tolerance
- Junction overlap ≥95%

**Small differences are OK:**
- Multi-mapper tie-breaking (order is arbitrary)
- MAPQ for multi-mappers (±5)
- Low-coverage junctions (1-2 read difference)

---

## Next Steps

1. **Run all tests** to establish baseline:
   ```bash
   ./run_tests.sh --all
   ```

2. **Save golden outputs** for passed tests:
   ```bash
   ./save_golden.sh yeast_100
   ./save_golden.sh yeast_1k
   # ... etc
   ```

3. **Integrate with workflow:**
   - Add `./ci.sh` to pre-commit hook
   - Run tests after code changes
   - Track alignment accuracy over time

4. **Expand coverage (future):**
   - Add larger test datasets
   - Add STARsolo test cases (Phase 14)
   - Add human genome tests

---

## Troubleshooting

### "STAR binary not found"
```bash
export STAR_BIN=/path/to/STAR
```

### "Test data not found"
Check that `test/data/small/yeast/` exists with reads and indices.

### "rustar-aligner binary not found"
```bash
cd .. && cargo build --release
```

### "Python script failed"
Check Python version:
```bash
python3 --version  # Should be 3.8+
```

---

## More Information

- **Full documentation:** `README.md`
- **STAR source reference:** `STAR_REFERENCE.md`
- **Implementation details:** `IMPLEMENTATION_SUMMARY.md`

---

**Questions?** Check the documentation or run `./verify_framework.sh` to diagnose issues.

**Ready to test?** Run `./run_tests.sh yeast_100` to get started!
