# PAL2NAL — Source Code Survey (2026-04-18)

## Background

PAL2NAL converts a multiple protein sequence alignment and the corresponding
DNA (or mRNA) sequences into a codon-based DNA alignment.  It was developed
by Mikita Suyama at the Bork Group (EMBL) and later at Kyoto University and
Kyushu University.

**Citation:** Mikita Suyama, David Torrents, and Peer Bork (2006).
PAL2NAL: robust conversion of protein sequence alignments into the
corresponding codon alignments.
*Nucleic Acids Res.* 34, W609–W612.
<https://www.bork.embl.de/pal2nal/pal2nal.pdf>

---

## Official website

| Item | Detail |
|---|---|
| Original URL | `http://www.bork.embl.de/pal2nal/` |
| Current status | **403 Forbidden** (as of 2026-04-18) |
| Wayback Machine captures | **147 captures** (earliest ~2006, latest 2025-11-19) |
| Latest archived snapshot | <https://web.archive.org/web/20251119055114/http://www.bork.embl.de/pal2nal/> |
| Download link on site | `http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz` |
| Download via Wayback | <https://web.archive.org/web/20251119055114/http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz> |

The "What's new" section on the official page **stops at December 2, 2011**
(v14).  No newer entry has ever appeared there.  The download tarball has been
byte-for-bit identical across all archived snapshots:

```
SHA-256: 0e458d38571ab6d61078047a9e80d843bb2e733414eeff8b1144d1c6bc13848f
File:    pal2nal.v14.tar.gz  (17031 bytes, internal mtime 2011-12-02)
```

---

## Version history (from official "What's new")

| Date | Version | Change |
|---|---|---|
| 2007-02-05 | v12 | Same-ID cross-referencing; speed improvements; `bl2seq` error helper |
| 2009-04-03 | — | Mirror added at Kyoto University |
| 2009-09-08 | — | S/N site counts added to dS/dN output |
| 2010-07-26 | v13 | Script update |
| **2011-12-02** | **v14** | All NCBI codon tables added (23 tables) — **last release on official site** |
| **2017-10-17** | **v14.1** | GPL-2.0 license header added — **never posted to official site** |

---

## Known public copies

### 1. Official tarball (v14, 2011)

- **Source:** `http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz`
  (accessible via Wayback Machine only)
- **Contents:** `pal2nal.pl`, `test.aln`, `test.nuc`, `README`,
  `for_paml/` (test.cnt, test.tree, test.codeml.ori)
- **Version string in script:** `pal2nal.pl  (v14)`
- **License:** none stated in the script itself

### 2. GitHub — liaochenlanruo/PAL2NAL (v14.1, 2017) ✅ **best available copy**

- **URL:** <https://github.com/liaochenlanruo/PAL2NAL>
- **Stars / forks:** 9 ⭐ / 2
- **Created:** 2019-06-02 (pushed 2019-06-02)
- **Provenance:** Mikita Suyama sent this version **directly by email** to
  liaochenlanruo for the purpose of writing a bioconda recipe.
- **Contents:** same as official tarball plus `LICENSE` (GPL-2.0) and
  `.gitignore`
- **Version string in script:** `pal2nal.pl  (v14.1)`
- **License:** GPL-2.0 (explicit in both `LICENSE` file and script header)
- **This is the copy cloned into this directory.**

### 3. GitHub — aur-archive/pal2nal

- **URL:** <https://github.com/aur-archive/pal2nal>
- **Contents:** Arch Linux AUR `PKGBUILD` only — **no script**
- **Stars:** 1

### 4. GitHub — dukecomeback/pal4nal.pl

- **URL:** <https://github.com/dukecomeback/pal4nal.pl>
- **Description:** "similar as what pal2nal.pl do … except some differences"
- **Note:** This is a **derivative work**, not a copy of the original

### 5. GitHub — josephryan/pal2nal_gblocker

- **URL:** <https://github.com/josephryan/pal2nal_gblocker>
- **Description:** Produces a nucleotide alignment matching a Gblocked
  protein alignment — a **downstream wrapper**, not the script itself

### 6. GitHub — nvk747/docker_pal2nal (and _CWL, _snakemake, _nextflow)

- Docker/CWL/Snakemake/Nextflow workflow wrappers — no standalone script

### 7. GitHub — inutano/pal2nal-cwl

- CWL wrapper (Dockerfile-based) — no standalone script

### 8. BIRCH bioinformatics system (Univ. of Manitoba)

- **Package page:** <https://home.cc.umanitoba.ca/~psgendb/birchdocdb/package/pal2nal.html>
- **Status:** BIRCH lists `pal2nal.pl` in its package index and points to
  `http://www.bork.embl.de/pal2nal/#Ref` for documentation.
  BIRCH bundles the upstream tarball unchanged — it has **no independent
  git repository** for pal2nal and does **no independent development**.
- No `pal2nal.py` exists anywhere (BIRCH or otherwise).

---

## Differences between v14 (official tarball) and v14.1 (this repo)

The diff is **purely cosmetic** — only comment lines changed, zero functional
difference:

```diff
--- pal2nal.v14/pal2nal.pl   (official tarball, 2011-12-02)
+++ PAL2NAL/pal2nal.pl       (liaochenlanruo GitHub, v14.1, 2017-10-17)

1a2
> #    pal2nal.pl (c) 2017 Mikita Suyama <mikita@bioreg.kyushu-u.ac.jp>

3c4,10
< #    pal2nal.pl  (v14)                                      Mikita Suyama
---
> #    #----------------------------------------#
> #      This script is licensed under GPL v2.0
> #    #----------------------------------------#
> #
> #
> #
> #    pal2nal.pl  (v14.1)                                      Mikita Suyama

58a66,69
> #
> #  v14.1
> #     - licensed under GPL v2.0
> #                                                          2017/10/17
```

**Summary of changes:**
- Copyright + contact email line added at top
- Explicit GPL v2.0 license notice block added
- `v14.1` changelog entry with date `2017/10/17` added at end of header

No logic, no algorithm, no codon table, no output format changed.

---

## Conclusion

**v14.1 in this repository is the definitive final release of PAL2NAL.**

There is no newer public version.  The only outstanding known issue is that
the official website (bork.embl.de) was never updated to advertise v14.1 or
to serve the updated tarball — the 2017 GPL-licensing revision exists only in
this GitHub repository.

The tool is **stable legacy code** with no active development.  For production
use, prefer this repository's `pal2nal.pl` (v14.1) over the official tarball
since it carries an explicit GPL-2.0 license, which is required for packaging
(bioconda, AUR, etc.).

---

## Contact

Mikita Suyama  
`mikita@bioreg.kyushu-u.ac.jp`  
(as of November 2025 Wayback snapshot)
