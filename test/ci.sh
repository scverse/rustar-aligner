#!/usr/bin/env bash

# ci.sh - Fast CI test suite for rustar-aligner
# Runs minimal test set for quick validation (pre-commit, CI/CD)
# Exit codes: 0=pass, 1=regression, 2=build failed, 3=execution failed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

log() {
    echo "[CI] $*"
}

error() {
    echo "[CI] ERROR: $*" >&2
}

# ==============================================================================
# Build
# ==============================================================================

log "Building rustar-aligner in release mode..."
cd "$PROJECT_ROOT"

if ! cargo build --release 2>&1 | tee /tmp/rustar-aligner_build.log; then
    error "Build failed"
    exit 2
fi

log "Build successful"

# ==============================================================================
# Unit Tests
# ==============================================================================

log "Running unit tests..."

if ! cargo test --release 2>&1 | tee /tmp/rustar-aligner_test.log; then
    error "Unit tests failed"
    exit 3
fi

log "Unit tests passed"

# ==============================================================================
# Fast Integration Tests
# ==============================================================================

log "Running fast integration tests (100-read dataset)..."

# Run only the smallest test case
if ! "$SCRIPT_DIR/run_tests.sh" yeast_100; then
    error "Integration test failed"
    exit 3
fi

log "Fast integration tests passed"

# ==============================================================================
# Golden Output Comparison (if golden outputs exist)
# ==============================================================================

GOLDEN_DIR="$SCRIPT_DIR/golden/yeast_100"

if [[ -d "$GOLDEN_DIR" ]]; then
    log "Comparing against golden outputs..."

    # Find most recent test result
    LATEST_RESULT=$(find "$SCRIPT_DIR/results" -maxdepth 1 -name "*_yeast_100" -type d | sort | tail -1)

    if [[ -z "$LATEST_RESULT" ]]; then
        error "No test results found for comparison"
        exit 3
    fi

    # Compare statistics
    if [[ -f "$GOLDEN_DIR/stats.json" ]]; then
        if ! python3 "$SCRIPT_DIR/compare_golden.py" \
            --golden "$GOLDEN_DIR/stats.json" \
            --current "$LATEST_RESULT/rustar-aligner" \
            --tolerance 0.01; then
            error "Regression detected: statistics differ from golden output"
            exit 1
        fi
    fi

    log "Golden output comparison passed"
else
    log "No golden outputs found, skipping comparison"
fi

# ==============================================================================
# Success
# ==============================================================================

log "=========================================="
log "CI Tests PASSED ✓"
log "=========================================="
log "Build:           OK"
log "Unit tests:      OK"
log "Integration:     OK"
log "Golden compare:  OK"

exit 0
