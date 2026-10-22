#!/usr/bin/env bash
# Generate THIRD_PARTY_LICENSES.md: the licence notices of every crate linked
# into rustar-aligner, which the release archives and the Docker image ship
# (their licences require the notices in binary copies too).
#
# Usage: scripts/third-party-licenses.sh [OUTPUT]   (default THIRD_PARTY_LICENSES.md)
#
# Downloads a pinned cargo-about release, checks its SHA-256, and runs it with
# the repo's about.toml / about.hbs. Linux x86_64 and aarch64 only (the CI,
# release and Docker builders). Fails if a dependency's licence is not in
# about.toml's accepted list, or a clarified licence file has changed.
set -euo pipefail

VERSION="0.9.2"
OUTPUT="${1:-THIRD_PARTY_LICENSES.md}"

case "$(uname -m)" in
    x86_64)
        ARCH="x86_64"
        SHA256="9099a59e820c38a68b9d65f300662a567d56562f9a10f6aa4c7e86c17c2566af"
        ;;
    aarch64 | arm64)
        ARCH="aarch64"
        SHA256="af5169282fb6f84e13471493f405437e43ac517744c9ae12fbe2cdf0a6f0e5a8"
        ;;
    *)
        echo "third-party-licenses: unsupported architecture $(uname -m)" >&2
        exit 1
        ;;
esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -r "$TMP"' EXIT

NAME="cargo-about-${VERSION}-${ARCH}-unknown-linux-musl"
curl -sSfL -o "$TMP/$NAME.tar.gz" \
    "https://github.com/EmbarkStudios/cargo-about/releases/download/${VERSION}/${NAME}.tar.gz"
echo "${SHA256}  $TMP/$NAME.tar.gz" | sha256sum -c --quiet -
tar xzf "$TMP/$NAME.tar.gz" -C "$TMP"

cd "$ROOT"
# cargo-about only warns when a clarified licence file no longer matches its
# checksum, and then drops that clarification (and with it the notice) from
# the output. A clean run logs nothing, so treat any warning as a failure.
"$TMP/$NAME/cargo-about" generate --locked --fail about.hbs -o "$OUTPUT" 2> "$TMP/log" || {
    cat "$TMP/log" >&2
    exit 1
}
if grep -qiE "warn|error" "$TMP/log"; then
    cat "$TMP/log" >&2
    echo "third-party-licenses: cargo-about warned; refusing an incomplete notice file" >&2
    rm -f "$OUTPUT"
    exit 1
fi
echo "third-party-licenses: wrote $OUTPUT"
