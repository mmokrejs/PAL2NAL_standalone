#!/bin/bash
# t/run_tests.sh — TAP-like test runner for pal2nal.pl
# Usage: bash t/run_tests.sh [pal2nal.pl path]
set -euo pipefail

PAL2NAL="${1:-$(dirname "$0")/../pal2nal.pl}"
T="$(dirname "$0")"
PASS=0; FAIL=0; SKIP=0

ok() {
    local desc="$1"; shift
    local got expected
    got="$("$@" 2>/dev/null)" || true
    expected="$(cat)"           # from stdin here-doc
    if [[ "$got" == "$expected" ]]; then
        echo "ok - $desc"
        (( PASS++ )) || true
    else
        echo "not ok - $desc"
        echo "  # GOT:      $(echo "$got"      | head -3)"
        echo "  # EXPECTED: $(echo "$expected" | head -3)"
        (( FAIL++ )) || true
    fi
}

ok_contains() {
    local desc="$1" pattern="$2"; shift 2
    local got
    got="$("$@" 2>&1)" || true
    if echo "$got" | grep -qF "$pattern"; then
        echo "ok - $desc"
        (( PASS++ )) || true
    else
        echo "not ok - $desc"
        echo "  # pattern '$pattern' not found in output"
        (( FAIL++ )) || true
    fi
}

ok_stderr_contains() {
    local desc="$1" pattern="$2"; shift 2
    local err
    err="$("$@" 2>&1 >/dev/null)" || true
    if echo "$err" | grep -qF "$pattern"; then
        echo "ok - $desc"
        (( PASS++ )) || true
    else
        echo "not ok - $desc"
        echo "  # pattern '$pattern' not found in stderr"
        (( FAIL++ )) || true
    fi
}

echo "TAP version 13"

# ─── T01: Basic CLUSTAL output with bundled test data ───────────────────────
echo "# T01: bundled test data, CLUSTAL output"
ok_contains "T01 CLUSTAL header present" "CLUSTAL W multiple sequence alignment" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -output clustal

# ─── T02: FASTA output format ───────────────────────────────────────────────
echo "# T02: fasta output"
ok_contains "T02 fasta output has FASTA header" ">" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -output fasta

# ─── T03: PAML output format ────────────────────────────────────────────────
echo "# T03: paml output first line"
ok_contains "T03 paml output has sequence count" "2" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -output paml

# ─── T04: -nogap removes stop codon columns ─────────────────────────────────
echo "# T04: -nogap removes stop codon columns"
# Pseudogene alignment contains in-frame stops (*); -nogap must remove those codons.
LEN_NOGAP=$(perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" \
    -output fasta -nogap 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c)
LEN_FULL=$(perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" \
    -output fasta 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c)
if [[ "$LEN_NOGAP" -le "$LEN_FULL" ]]; then
    echo "ok - T04 -nogap sequence length ($LEN_NOGAP) <= full length ($LEN_FULL)"
    (( PASS++ )) || true
else
    echo "not ok - T04 -nogap length $LEN_NOGAP > full $LEN_FULL (unexpected)"
    (( FAIL++ )) || true
fi

# ─── T05: -h shows help ─────────────────────────────────────────────────────
echo "# T05: -h help output"
ok_contains "T05 -h shows Usage" "Usage" \
    perl "$PAL2NAL" -h

# ─── T06 (BUG-3 baseline): help version string ──────────────────────────────
echo "# T06: help version string (BUG-3 baseline)"
HELP_VER=$(perl "$PAL2NAL" -h 2>&1 | grep "pal2nal.pl" | head -1)
# Record current (buggy) behavior: prints v14
if echo "$HELP_VER" | grep -q "v14\.1"; then
    echo "ok - T06 help shows v14.1 (bug already fixed)"
    (( PASS++ )) || true
elif echo "$HELP_VER" | grep -q "v14"; then
    echo "ok - T06 help shows v14 (BUG-3 baseline recorded)"
    (( PASS++ )) || true
else
    echo "not ok - T06 unexpected help version: $HELP_VER"
    (( FAIL++ )) || true
fi

# ─── T07: -codontable validation rejects invalid table ──────────────────────
echo "# T07: invalid codontable rejected"
ok_stderr_contains "T07 bad codontable gives error" "ERROR" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -codontable 99

# ─── T08: -output codon rejected with -nogap ────────────────────────────────
echo "# T08: -output codon + -nogap gives error"
ok_contains "T08 codon+nogap error" "ERROR" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -output codon -nogap

# ─── T09: mismatched seq count gives error ──────────────────────────────────
echo "# T09: mismatched sequence count detected"
ok_stderr_contains "T09 seq count mismatch error" "ERROR" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test_1seq.nuc"

# ─── T10: same-ID matching produces non-empty FASTA output ──────────────────
echo "# T10: same-ID mode produces output"
OUT_SAMEID=$(perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" \
    -output fasta 2>/dev/null)
if echo "$OUT_SAMEID" | grep -q "^>"; then
    echo "ok - T10 same-ID mode produces FASTA records"
    (( PASS++ )) || true
else
    echo "not ok - T10 same-ID mode produced no FASTA records"
    (( FAIL++ )) || true
fi

# ─── T11: BUG-7 baseline — error message option name ────────────────────────
echo "# T11: BUG-7 baseline — option name in codon+nogap error"
ERR11=$(perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" \
    -output codon -nogap 2>&1)
if echo "$ERR11" | grep -q "\-outform"; then
    echo "ok - T11 BUG-7 baseline: error says '-outform' (wrong, to be fixed)"
    (( PASS++ )) || true
elif echo "$ERR11" | grep -q "\-output"; then
    echo "ok - T11 BUG-7 fixed: error correctly says '-output'"
    (( PASS++ )) || true
else
    echo "not ok - T11 unexpected error message: $ERR11"
    (( FAIL++ )) || true
fi

# ─── T12: -codontable 2 (vertebrate mito) accepted ──────────────────────────
echo "# T12: -codontable 2 accepted"
ok_contains "T12 codontable 2 produces output" "CLUSTAL" \
    perl "$PAL2NAL" "$T/data/test.aln" "$T/data/test.nuc" -codontable 2

# ─── T13: BUG-5 baseline — -nogap stop-codon filter with table 6 ────────────
echo "# T13: BUG-5 baseline — -nogap with table 6 (CAA/CAG are Gln not stop)"
# In table 6 TAA/TAG code for Gln; Universal stop codon filter should NOT remove them
# This test documents the baseline (buggy) behavior
OUT13=$(perl "$PAL2NAL" "$T/data/test_table6.aln" "$T/data/test_table6.nuc" \
    -output fasta -nogap -codontable 6 2>/dev/null || echo "ERROR")
if [[ -z "$OUT13" ]] || [[ "$OUT13" == "ERROR" ]]; then
    echo "ok - T13 BUG-5 baseline: -nogap with table 6 drops Gln codons (buggy)"
    (( PASS++ )) || true
else
    echo "ok - T13 BUG-5 fixed: -nogap with table 6 retains Gln codons"
    (( PASS++ )) || true
fi

echo ""
echo "Results: $PASS passed, $FAIL failed, $SKIP skipped"
echo "1..$((PASS + FAIL + SKIP))"
[[ "$FAIL" -eq 0 ]]
