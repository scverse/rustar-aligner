#!/usr/bin/env bash

# run_tests.sh - Master test orchestration script for rustar-aligner vs STAR comparison
# Usage: ./run_tests.sh [--all | test_name1 test_name2 ...] [--parallel] [--keep-all]

set -euo pipefail

# ==============================================================================
# Configuration
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RUSTAR_BIN="$PROJECT_ROOT/target/release/rustar-aligner-aligner"
STAR_BIN="${STAR_BIN:-$(which STAR || echo "")}"

# Test data paths
DATA_DIR="$SCRIPT_DIR/data/small/yeast"
READS_DIR="$DATA_DIR/reads"
STAR_GENOME_DIR="$DATA_DIR/indices"
RUSTAR_ALIGNER_GENOME_DIR="$DATA_DIR/indices_rustar"
GTF_FILE="$DATA_DIR/reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf"

# Output paths
RESULTS_DIR="$SCRIPT_DIR/results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Default parameters
PARALLEL=false
KEEP_ALL=false
TIMEOUT=600  # 10 minutes per test
SELECTED_TESTS=()

# ==============================================================================
# Test Case Definitions
# ==============================================================================

# Format: "name:dataset:reads:mode:extra_args"
# mode: single|paired
TEST_CASES=(
    "yeast_100:yeast:ERR12389696_sub_1_100.fastq.gz:single:--outSAMtype SAM"
    "yeast_1k:yeast:ERR12389696_sub_1_1000.fastq.gz:single:--outSAMtype SAM"
    "yeast_1k_bam:yeast:ERR12389696_sub_1_1000.fastq.gz:single:--outSAMtype BAM Unsorted"
    "yeast_1k_paired:yeast:ERR12389696_sub_1_1000.fastq.gz,ERR12389696_sub_2_1000.fastq.gz:paired:--outSAMtype SAM"
    "yeast_1k_twopass:yeast:ERR12389696_sub_1_1000.fastq.gz:single:--twopassMode Basic --outSAMtype SAM"
    "yeast_10k:yeast:ERR12389696_sub_1_10k.fastq.gz:single:--outSAMtype SAM"
)

# ==============================================================================
# Helper Functions
# ==============================================================================

print_usage() {
    cat << EOF
Usage: $0 [OPTIONS] [TEST_NAMES...]

Options:
    --all           Run all test cases
    --parallel      Run test cases in parallel
    --keep-all      Keep all intermediate files
    --timeout SEC   Timeout per test (default: 600)
    --help          Show this help

Test cases:
EOF
    for test in "${TEST_CASES[@]}"; do
        local name="${test%%:*}"
        echo "    $name"
    done
    echo ""
    echo "Examples:"
    echo "    $0 yeast_100              # Run single test"
    echo "    $0 yeast_100 yeast_1k     # Run multiple tests"
    echo "    $0 --all                  # Run all tests"
    echo "    $0 --all --parallel       # Run all tests in parallel"
}

log() {
    echo "[$(date +%T)] $*"
}

error() {
    echo "[$(date +%T)] ERROR: $*" >&2
}

check_prerequisites() {
    log "Checking prerequisites..."

    # Check rustar-aligner binary
    if [[ ! -x "$RUSTAR_BIN" ]]; then
        error "rustar-aligner binary not found at $RUSTAR_BIN"
        error "Run: cargo build --release"
        exit 1
    fi

    # Check STAR binary
    if [[ -z "$STAR_BIN" ]]; then
        error "STAR binary not found in PATH"
        error "Install STAR or set STAR_BIN environment variable"
        exit 1
    fi

    # Check Python and dependencies
    if ! command -v python3 &> /dev/null; then
        error "Python 3 not found"
        exit 1
    fi

    # Check test data
    if [[ ! -d "$DATA_DIR" ]]; then
        error "Test data not found at $DATA_DIR"
        error "Run: cd test/data/small/yeast && ./test_yeast.sh setup"
        exit 1
    fi

    log "Prerequisites OK"
}

# ==============================================================================
# Test Execution
# ==============================================================================

run_star() {
    local output_dir="$1"
    local reads="$2"
    local extra_args="$3"

    mkdir -p "$output_dir"

    log "Running STAR..."

    # Parse extra args to extract relevant flags
    local twopass_mode=""
    if [[ "$extra_args" == *"--twopassMode Basic"* ]]; then
        twopass_mode="--twopassMode Basic"
    fi

    local sam_type="SAM"
    if [[ "$extra_args" == *"BAM"* ]]; then
        sam_type="BAM Unsorted"
    fi

    # Build STAR command
    local cmd=(
        "$STAR_BIN"
        --runMode alignReads
        --genomeDir "$STAR_GENOME_DIR"
        --readFilesIn ${reads//,/ }
        --readFilesCommand zcat
        --outFileNamePrefix "$output_dir/"
        --outSAMtype $sam_type
        --sjdbGTFfile "$GTF_FILE"
        --runThreadN 4
    )

    if [[ -n "$twopass_mode" ]]; then
        cmd+=(--twopassMode Basic)
    fi

    # Run with timeout
    if timeout "$TIMEOUT" "${cmd[@]}" > "$output_dir/star.log" 2>&1; then
        log "STAR completed successfully"
        return 0
    else
        error "STAR failed or timed out"
        cat "$output_dir/star.log" >&2
        return 1
    fi
}

run_rustar_aligner() {
    local output_dir="$1"
    local reads="$2"
    local extra_args="$3"

    mkdir -p "$output_dir"

    log "Running rustar-aligner..."

    # Build rustar-aligner command
    local cmd=(
        "$RUSTAR_BIN"
        --runMode alignReads
        --genomeDir "$RUSTAR_ALIGNER_GENOME_DIR"
        --readFilesIn ${reads//,/ }
        --readFilesCommand zcat
        --outFileNamePrefix "$output_dir/"
        --sjdbGTFfile "$GTF_FILE"
        --runThreadN 4
    )

    # Add extra args
    if [[ -n "$extra_args" ]]; then
        eval "cmd+=($extra_args)"
    fi

    # Run with timeout
    if timeout "$TIMEOUT" "${cmd[@]}" > "$output_dir/rustar-aligner-aligner.log" 2>&1; then
        log "rustar-aligner completed successfully"
        return 0
    else
        error "rustar-aligner failed or timed out"
        cat "$output_dir/rustar-aligner-aligner.log" >&2
        return 1
    fi
}

compare_outputs() {
    local test_dir="$1"
    local star_dir="$test_dir/star"
    local rustar_aligner_dir="$test_dir/rustar-aligner-aligner"
    local comparison_dir="$test_dir/comparison"

    mkdir -p "$comparison_dir"

    log "Comparing outputs..."

    local status=0

    # Compare SAM/BAM files
    if [[ -f "$star_dir/Aligned.out.sam" && -f "$rustar_aligner_dir/Aligned.out.sam" ]]; then
        if python3 "$SCRIPT_DIR/compare_sam.py" \
            --star "$star_dir/Aligned.out.sam" \
            --rustar-aligner "$rustar_aligner_dir/Aligned.out.sam" \
            --tolerance 0.01 \
            --output "$comparison_dir/alignment_diff.txt" \
            > "$comparison_dir/sam_comparison.log" 2>&1; then
            log "SAM comparison: PASS"
        else
            error "SAM comparison: FAIL"
            status=1
        fi
    elif [[ -f "$star_dir/Aligned.out.bam" && -f "$rustar_aligner_dir/Aligned.out.bam" ]]; then
        # For BAM, convert to SAM first
        log "Converting BAM to SAM for comparison..."
        samtools view -h "$star_dir/Aligned.out.bam" > "$comparison_dir/star_tmp.sam"
        samtools view -h "$rustar_aligner_dir/Aligned.out.bam" > "$comparison_dir/rustar_aligner_tmp.sam"

        if python3 "$SCRIPT_DIR/compare_sam.py" \
            --star "$comparison_dir/star_tmp.sam" \
            --rustar-aligner "$comparison_dir/rustar_aligner_tmp.sam" \
            --tolerance 0.01 \
            --output "$comparison_dir/alignment_diff.txt" \
            > "$comparison_dir/sam_comparison.log" 2>&1; then
            log "BAM comparison: PASS"
        else
            error "BAM comparison: FAIL"
            status=1
        fi

        rm -f "$comparison_dir/star_tmp.sam" "$comparison_dir/rustar_aligner_tmp.sam"
    fi

    # Compare junction files
    if [[ -f "$star_dir/SJ.out.tab" && -f "$rustar_aligner_dir/SJ.out.tab" ]]; then
        if python3 "$SCRIPT_DIR/compare_junctions.py" \
            --star "$star_dir/SJ.out.tab" \
            --rustar-aligner "$rustar_aligner_dir/SJ.out.tab" \
            --tolerance 0.10 \
            --output "$comparison_dir/junction_diff.txt" \
            > "$comparison_dir/junction_comparison.log" 2>&1; then
            log "Junction comparison: PASS"
        else
            error "Junction comparison: FAIL"
            status=1
        fi
    fi

    # Generate summary
    {
        echo "=========================================="
        echo "Test Comparison Summary"
        echo "=========================================="
        echo "Test directory: $test_dir"
        echo "Timestamp: $(date)"
        echo ""

        if [[ -f "$comparison_dir/alignment_diff.txt" ]]; then
            echo "--- Alignment Comparison ---"
            cat "$comparison_dir/alignment_diff.txt"
            echo ""
        fi

        if [[ -f "$comparison_dir/junction_diff.txt" ]]; then
            echo "--- Junction Comparison ---"
            cat "$comparison_dir/junction_diff.txt"
            echo ""
        fi

        if [[ $status -eq 0 ]]; then
            echo "Overall Status: PASS ✓"
        else
            echo "Overall Status: FAIL ✗"
        fi
    } > "$comparison_dir/summary.txt"

    cat "$comparison_dir/summary.txt"

    return $status
}

run_test_case() {
    local test_spec="$1"

    # Parse test specification
    IFS=':' read -r name dataset reads mode extra_args <<< "$test_spec"

    local test_dir="$RESULTS_DIR/${TIMESTAMP}_${name}"
    local star_dir="$test_dir/star"
    local rustar_aligner_dir="$test_dir/rustar-aligner-aligner"

    log "=========================================="
    log "Running test: $name"
    log "=========================================="

    # Prepare read paths
    local read_paths=""
    IFS=',' read -ra READ_FILES <<< "$reads"
    for read_file in "${READ_FILES[@]}"; do
        if [[ -n "$read_paths" ]]; then
            read_paths="$read_paths,$READS_DIR/$read_file"
        else
            read_paths="$READS_DIR/$read_file"
        fi
    done

    # Check if reads exist
    IFS=',' read -ra READ_PATHS <<< "$read_paths"
    for read_path in "${READ_PATHS[@]}"; do
        if [[ ! -f "$read_path" ]]; then
            error "Read file not found: $read_path"
            return 1
        fi
    done

    # Run STAR
    if ! run_star "$star_dir" "$read_paths" "$extra_args"; then
        error "STAR execution failed for $name"
        echo "FAILED" > "$test_dir/FAILED"
        return 1
    fi

    # Run rustar-aligner
    if ! run_rustar_aligner "$rustar_aligner_dir" "$read_paths" "$extra_args"; then
        error "rustar-aligner execution failed for $name"
        echo "FAILED" > "$test_dir/FAILED"
        return 1
    fi

    # Compare outputs
    if compare_outputs "$test_dir"; then
        log "Test $name: PASS ✓"
        echo "PASSED" > "$test_dir/PASSED"

        # Cleanup intermediate files if not keeping all
        if [[ "$KEEP_ALL" == false ]]; then
            rm -f "$star_dir"/_STARtmp "$star_dir"/Log.* "$star_dir"/SJ.out.tab
            rm -f "$rustar_aligner_dir"/Log.*
        fi

        return 0
    else
        error "Test $name: FAIL ✗"
        echo "FAILED" > "$test_dir/FAILED"
        return 1
    fi
}

# ==============================================================================
# Main Execution
# ==============================================================================

main() {
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --help|-h)
                print_usage
                exit 0
                ;;
            --all)
                SELECTED_TESTS=("${TEST_CASES[@]}")
                shift
                ;;
            --parallel)
                PARALLEL=true
                shift
                ;;
            --keep-all)
                KEEP_ALL=true
                shift
                ;;
            --timeout)
                TIMEOUT="$2"
                shift 2
                ;;
            -*)
                error "Unknown option: $1"
                print_usage
                exit 1
                ;;
            *)
                # Find matching test case
                local found=false
                for test in "${TEST_CASES[@]}"; do
                    if [[ "${test%%:*}" == "$1" ]]; then
                        SELECTED_TESTS+=("$test")
                        found=true
                        break
                    fi
                done
                if [[ "$found" == false ]]; then
                    error "Unknown test case: $1"
                    exit 1
                fi
                shift
                ;;
        esac
    done

    # Default to all tests if none selected
    if [[ ${#SELECTED_TESTS[@]} -eq 0 ]]; then
        SELECTED_TESTS=("${TEST_CASES[@]}")
    fi

    check_prerequisites

    mkdir -p "$RESULTS_DIR"

    log "Starting test suite with ${#SELECTED_TESTS[@]} test(s)"

    local passed=0
    local failed=0
    local pids=()

    # Run tests
    for test in "${SELECTED_TESTS[@]}"; do
        if [[ "$PARALLEL" == true ]]; then
            run_test_case "$test" &
            pids+=($!)
        else
            if run_test_case "$test"; then
                ((passed++))
            else
                ((failed++))
            fi
        fi
    done

    # Wait for parallel tests
    if [[ "$PARALLEL" == true ]]; then
        for pid in "${pids[@]}"; do
            if wait "$pid"; then
                ((passed++))
            else
                ((failed++))
            fi
        done
    fi

    # Print summary
    log "=========================================="
    log "Test Suite Complete"
    log "=========================================="
    log "Passed: $passed"
    log "Failed: $failed"
    log "Total:  $((passed + failed))"
    log "Results directory: $RESULTS_DIR"

    if [[ $failed -eq 0 ]]; then
        log "All tests PASSED ✓"
        exit 0
    else
        error "Some tests FAILED ✗"
        exit 1
    fi
}

main "$@"
