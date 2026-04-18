# PAL2NAL v14.1 — Bug Report & Fixes (2026-04-18)

Full code read of `pal2nal.pl` (2010 lines). Issues are ranked by severity and have been resolved via individual Pull Requests against the mirrored repository.

---

## 10 bugs found across 4 severity levels

| ID | Severity | Bug | PR |
|---|---|---|---|
| **BUG-1** | 🔴 High | **Wrong variable in frameshift column sizing** — the `-output codon` path computes `$tmplen` from `$tmpaa` (last-sequence value) instead of `$maxaan` (the maximum across all sequences). Produces wrong codon column widths when sequences have different frameshift lengths at the same position. | [#7](https://github.com/mmokrejs/PAL2NAL/pull/7) |
| **BUG-2** | 🔴 High | **`peppos` counter wrong in anchor fallback loop** — `$i * 10 + $j` doesn't account for gap characters inside anchor chunks, so mismatch warning messages report wrong alignment positions. *(Note: only affects a commented-out debug line).* | [#8](https://github.com/mmokrejs/PAL2NAL/pull/8) |
| **BUG-3** | 🟡 Medium | **Help text says `v14`**, not `v14.1` — missed when the license-only bump was applied in 2017. | [#1](https://github.com/mmokrejs/PAL2NAL/pull/1) |
| **BUG-4** | 🟠 Medium | **Shell injection + temp-file race** in the `bl2seq` diagnostic path — sequence IDs are interpolated directly into `system()` unescaped. Replaced with safe instruction to use modern `blastp`. | [#9](https://github.com/mmokrejs/PAL2NAL/pull/9) |
| **BUG-5** | 🟡 Medium | **`-nogap` stop-codon regex is hardcoded Universal** (`TAA/TAG/TGA`) regardless of `-codontable`, so with e.g. table 6 it wrongly filters Gln codons. Fixed to respect `-codontable`. | [#6](https://github.com/mmokrejs/PAL2NAL/pull/6) |
| **BUG-6** | 🟡 Medium | **Gblocks format parser is broken** — `$getaln` is reset to `0` on every line before the data-processing branch is tested, so alignment data in Gblocks format is silently dropped. | [#5](https://github.com/mmokrejs/PAL2NAL/pull/5) |
| **BUG-7** | 🟢 Low | Error message says `-outform` but the real flag is `-output`. | [#2](https://github.com/mmokrejs/PAL2NAL/pull/2) |
| **BUG-8** | 🟢 Low | Typo: "alignemt" in comment. | [#3](https://github.com/mmokrejs/PAL2NAL/pull/3) |
| **BUG-9** | 🟢 Low | Typo: "exlucded" in commented-out line. | [#4](https://github.com/mmokrejs/PAL2NAL/pull/4) |
| **BUG-10** | 🟢 Low | No `use strict; use warnings;` — deprecated bareword filehandles and typeglob passing throughout. *(Unfixed: Requires major rewrite of legacy Perl 4 idioms)* | N/A |
