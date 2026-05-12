#!/usr/bin/env bash

# verify_framework.sh - Quick verification of test framework setup
# Checks prerequisites and validates file structure

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

pass() {
    echo -e "${GREEN}✓${NC} $*"
}

fail() {
    echo -e "${RED}✗${NC} $*"
}

warn() {
    echo -e "${YELLOW}⚠${NC} $*"
}

check_count=0
pass_count=0
fail_count=0

check() {
    ((check_count++))
    if "$@"; then
        ((pass_count++))
        return 0
    else
        ((fail_count++))
        return 1
    fi
}

echo "=========================================="
echo "rustar-aligner Test Framework Verification"
echo "=========================================="
echo ""

# Check executable scripts
echo "Checking test scripts..."
for script in run_tests.sh ci.sh investigate.sh save_golden.sh; do
    if check test -x "$SCRIPT_DIR/$script"; then
        pass "$script is executable"
    else
        fail "$script is NOT executable or missing"
    fi
done

# Check Python scripts
echo ""
echo "Checking comparison utilities..."
for script in compare_sam.py compare_junctions.py compare_chimeric.py compare_golden.py; do
    if check test -x "$SCRIPT_DIR/$script"; then
        pass "$script is executable"
    else
        fail "$script is NOT executable or missing"
    fi
done

# Check documentation
echo ""
echo "Checking documentation..."
for doc in README.md STAR_REFERENCE.md; do
    if check test -f "$SCRIPT_DIR/$doc"; then
        pass "$doc exists"
    else
        fail "$doc is missing"
    fi
done

# Check rustar-aligner binary
echo ""
echo "Checking rustar-aligner build..."
if check test -x "$PROJECT_ROOT/target/release/rustar-aligner"; then
    pass "rustar-aligner binary exists"
    version=$("$PROJECT_ROOT/target/release/rustar-aligner" --version 2>&1 || echo "unknown")
    echo "  Version: $version"
else
    fail "rustar-aligner binary NOT found"
    warn "Run: cargo build --release"
fi

# Check STAR
echo ""
echo "Checking STAR installation..."
if STAR_BIN=$(which STAR 2>/dev/null); then
    pass "STAR found at: $STAR_BIN"
    star_version=$(STAR --version 2>&1 | head -1 || echo "unknown")
    echo "  Version: $star_version"
else
    fail "STAR NOT found in PATH"
    warn "Install STAR or set STAR_BIN environment variable"
fi

# Check Python
echo ""
echo "Checking Python..."
if check command -v python3 &> /dev/null; then
    pass "Python 3 found"
    py_version=$(python3 --version)
    echo "  Version: $py_version"
else
    fail "Python 3 NOT found"
fi

# Check test data
echo ""
echo "Checking test data..."
DATA_DIR="$SCRIPT_DIR/data/small/yeast"
READS_DIR="$DATA_DIR/reads"

if check test -d "$DATA_DIR"; then
    pass "Test data directory exists"
else
    fail "Test data directory NOT found"
    warn "Run: cd test/data/small/yeast && ./test_yeast.sh setup"
fi

if check test -d "$READS_DIR"; then
    pass "Reads directory exists"

    # Count FASTQ files
    fastq_count=$(find "$READS_DIR" -name "*.fastq.gz" 2>/dev/null | wc -l)
    echo "  Found $fastq_count FASTQ files"

    # Check for specific test files
    for file in ERR12389696_sub_1_100.fastq.gz ERR12389696_sub_1_1000.fastq.gz ERR12389696_sub_1_10k.fastq.gz; do
        if test -f "$READS_DIR/$file"; then
            pass "  $file exists"
        else
            warn "  $file NOT found"
        fi
    done
else
    fail "Reads directory NOT found"
fi

if check test -d "$DATA_DIR/indices"; then
    pass "Genome index exists"
else
    fail "Genome index NOT found"
    warn "Run: cd test/data/small/yeast && ./test_yeast.sh setup"
fi

if check test -f "$DATA_DIR/reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf"; then
    pass "GTF annotation file exists"
else
    fail "GTF annotation file NOT found"
fi

# Check directory structure
echo ""
echo "Checking directory structure..."

if check test -d "$SCRIPT_DIR/golden"; then
    pass "Golden output directory exists"
else
    warn "Golden output directory NOT found (will be created on first save)"
fi

if test -d "$SCRIPT_DIR/results"; then
    warn "Results directory exists (contains old test runs)"
else
    pass "Results directory clean (will be created on first run)"
fi

if test -d "$SCRIPT_DIR/debug"; then
    warn "Debug directory exists (contains old debug files)"
else
    pass "Debug directory clean (will be created on first investigation)"
fi

# Summary
echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Total checks: $check_count"
pass "Passed: $pass_count"
if [[ $fail_count -gt 0 ]]; then
    fail "Failed: $fail_count"
else
    echo -e "${GREEN}Failed: 0${NC}"
fi

echo ""

if [[ $fail_count -eq 0 ]]; then
    pass "Test framework is ready!"
    echo ""
    echo "Try running:"
    echo "  ./run_tests.sh yeast_100          # Run single test"
    echo "  ./run_tests.sh --all              # Run all tests"
    echo "  ./ci.sh                           # Run fast CI tests"
    exit 0
else
    fail "Test framework has issues"
    echo ""
    echo "Please fix the issues above before running tests."
    exit 1
fi
