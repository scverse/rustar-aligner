//! Runtime CPU feature detection and compatibility guard for SIMD-tuned
//! binaries. The compile-time build flavour (OS, arch, SIMD level, git hash,
//! build timestamp) is assembled into `VERSION_BODY` by `build.rs` and read
//! via `env!("VERSION_BODY")`.

use anyhow::{Result, bail};

// ============================================================================
// Binary target (compile-time)
// ============================================================================

/// Returns the microarchitecture level this binary was compiled for.
pub fn binary_target() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    {
        if cfg!(target_feature = "avx512f") {
            "x86-64-v4"
        } else if cfg!(target_feature = "avx2") {
            "x86-64-v3"
        } else {
            "x86-64"
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if cfg!(target_feature = "sve") {
            "aarch64-neoverse-v1"
        } else {
            "aarch64"
        }
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        "unknown"
    }
}

// ============================================================================
// Runtime CPU feature detection
// ============================================================================

/// Detects CPU features available at runtime, ordered most → least capable.
pub fn detected_features() -> Vec<&'static str> {
    let mut features = Vec::new();

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512f") {
            features.push("AVX-512");
        }
        if is_x86_feature_detected!("avx2") {
            features.push("AVX2");
        }
        if is_x86_feature_detected!("sse4.2") {
            features.push("SSE4.2");
        }
        if is_x86_feature_detected!("popcnt") {
            features.push("POPCNT");
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("sve2") {
            features.push("SVE2");
        }
        if std::arch::is_aarch64_feature_detected!("sve") {
            features.push("SVE");
        }
        features.push("NEON");
    }

    features
}

// ============================================================================
// Upgrade hint
// ============================================================================

/// Returns an upgrade suggestion if a faster binary variant is available for
/// the current CPU, `None` if this binary is already optimal.
pub fn upgrade_hint() -> Option<String> {
    #[cfg(target_arch = "x86_64")]
    {
        let has_avx512 = is_x86_feature_detected!("avx512f");
        let has_avx2 = is_x86_feature_detected!("avx2");
        let target = binary_target();

        if target == "x86-64" && has_avx512 {
            return Some("Use the x86-64-v4 build for best performance.".to_string());
        }
        if target == "x86-64" && has_avx2 {
            return Some("Use the x86-64-v3 build for better performance.".to_string());
        }
        if target == "x86-64-v3" && has_avx512 {
            return Some("Use the x86-64-v4 build for best performance.".to_string());
        }
    }

    #[cfg(all(target_arch = "aarch64", target_os = "linux"))]
    {
        let has_sve = std::arch::is_aarch64_feature_detected!("sve");
        let target = binary_target();

        if target == "aarch64" && has_sve {
            return Some("Use the aarch64-neoverse-v1 build for better performance.".to_string());
        }
    }

    // No runtime hint for macOS aarch64: the `apple-m1` variant is a compiler
    // tuning target, not a distinct CPU feature set, so runtime detection
    // can't tell the two builds apart.

    None
}

// ============================================================================
// CPU compatibility guard
// ============================================================================

/// Checks that the CPU supports the features required by this binary.
/// Call at the very start of `main()` — before any SIMD code runs — to give
/// a friendly error instead of a SIGILL crash.
pub fn check_cpu_compat() -> Result<()> {
    #[cfg(target_arch = "x86_64")]
    match binary_target() {
        "x86-64-v4" if !is_x86_feature_detected!("avx512f") => bail!(
            "This rustar-aligner binary was compiled for x86-64-v4 (AVX-512) but your CPU \
             does not support AVX-512.\nPlease use the x86-64-v3 or baseline build instead.\n\
             See: https://github.com/Psy-Fer/rustar-aligner#installation"
        ),
        "x86-64-v3" if !is_x86_feature_detected!("avx2") => bail!(
            "This rustar-aligner binary was compiled for x86-64-v3 (AVX2) but your CPU \
             does not support AVX2.\nPlease use the baseline build instead.\n\
             See: https://github.com/Psy-Fer/rustar-aligner#installation"
        ),
        _ => {}
    }

    #[cfg(target_arch = "aarch64")]
    if binary_target() == "aarch64-neoverse-v1" && !std::arch::is_aarch64_feature_detected!("sve") {
        bail!(
            "This rustar-aligner binary was compiled for aarch64-neoverse-v1 (SVE) but your CPU \
             does not support SVE.\nPlease use the baseline aarch64 build instead.\n\
             See: https://github.com/Psy-Fer/rustar-aligner#installation"
        );
    }

    Ok(())
}

/// One-line summary of CPU features detected at runtime. No upgrade hint.
/// Example: `CPU: AVX2 SSE4.2 POPCNT detected`
pub fn cpu_detected_line() -> String {
    let features = detected_features();
    if features.is_empty() {
        "CPU: none detected".to_string()
    } else {
        format!("CPU: {} detected", features.join(" "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_target_is_not_empty() {
        assert!(!binary_target().is_empty());
    }

    #[test]
    fn version_body_contains_expected_fields() {
        let body = env!("VERSION_BODY");
        assert!(body.contains(env!("GIT_SHORT_HASH")));
        assert!(body.contains("built "));
    }

    #[test]
    fn cpu_detected_line_starts_with_cpu() {
        assert!(cpu_detected_line().starts_with("CPU:"));
    }

    #[test]
    fn check_cpu_compat_passes_on_current_hardware() {
        assert!(check_cpu_compat().is_ok());
    }
}
