# Deliberate divergences from STAR 2.7.11b

rustar-aligner aims to be a faithful port: given the same inputs and flags, it should
produce the same bytes as STAR 2.7.11b. This file records the places where it
deliberately does not, and why.

Every entry here is a case where STAR's behaviour is undefined, plainly wrong, or
impossible to reproduce, and where reproducing it would mean shipping a known
defect. Each is covered by a test that asserts the **correct** result, never
STAR's.

Divergences that arise from bugs on our side are not listed here. They are bugs,
and they get fixed.

---

## Format

Each entry gives: what STAR does, what rustar-aligner does instead, why, and the
test that locks the behaviour in.

---

## D-01 · Non-canonical annotated junctions carry no strand

**STAR.** Tracks a per-junction strand (`sjStr`) alongside the motif, and for an
annotated junction takes that strand from `sjdbStrand`. SJ.out.tab column 4 then
reports the annotated strand even when the motif is non-canonical.

**rustar-aligner.** Has no per-junction strand on the working transcript, so
column 4 reports `0` (undefined) for a non-canonical annotated junction where
STAR would report `1` or `2`. Canonical junctions are unaffected: their strand
is implied by the motif, which both derive identically.

**Why.** Not a deliberate improvement, just not yet ported. Recorded here so the
difference is not mistaken for a bug in the annotated-junction path, which is
otherwise faithful. Adding `junction_strands` to `WorkingTranscript` and
`Transcript` would close it.

**Test.** None yet; this entry is the marker.

---

## Divergences considered and rejected

Cases where STAR-rs, the other Rust port, made a different choice that
rustar-aligner deliberately does not adopt.

### Transcriptome-BAM primary flag

STAR picks the primary alignment for the transcriptome BAM with a per-thread
RNG, so the output depends on `--runThreadN`. STAR-rs resolves this by always
taking the first alignment (`j == 0`).

rustar-aligner already picks by a per-read seed, which is thread-count invariant
and therefore fixes the same defect. The two ports resolve it differently and
both diverge from STAR; there is no fidelity argument for switching to `j == 0`,
so rustar-aligner keeps its own rule.
