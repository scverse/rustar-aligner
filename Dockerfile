# ---- Build stage ----
FROM rust:1-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

ARG GIT_SHORT_HASH=unknown
ARG CPU_TARGET=""

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY src/ src/

RUN HOST_TRIPLE=$(rustc -vV | awk '/^host:/ {print $2}') && \
    GIT_SHORT_HASH="${GIT_SHORT_HASH}" \
    cargo build --release --target "$HOST_TRIPLE" \
        ${CPU_TARGET:+--config "target.'$HOST_TRIPLE'.rustflags=['-C', 'target-cpu=$CPU_TARGET']"} \
    && strip "target/$HOST_TRIPLE/release/rustar-aligner" \
    && cp "target/$HOST_TRIPLE/release/rustar-aligner" /rustar-aligner

# Licence notices of every crate linked into the binary (see about.toml).
COPY scripts/third-party-licenses.sh scripts/
COPY about.toml about.hbs LICENSE ./
COPY LICENSES/ LICENSES/
RUN scripts/third-party-licenses.sh /THIRD_PARTY_LICENSES.md

# ---- Runtime stage ----
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /rustar-aligner /usr/local/bin/rustar-aligner
COPY --from=builder /rustar-aligner /usr/local/bin/STAR
# rustar-aligner's licence and the third-party notices required in binary
# distributions: ported code (LICENSES/) and linked crates.
COPY LICENSE /usr/share/doc/rustar-aligner/LICENSE
COPY LICENSES/ /usr/share/doc/rustar-aligner/LICENSES/
COPY --from=builder /THIRD_PARTY_LICENSES.md /usr/share/doc/rustar-aligner/THIRD_PARTY_LICENSES.md

CMD ["rustar-aligner"]
