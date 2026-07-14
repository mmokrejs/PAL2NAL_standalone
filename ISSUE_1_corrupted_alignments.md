# Issue: pal2nal produces highly divergent/corrupted alignments on complex frameshifted sequences

## Description
When running `pal2nal.pl` on sequences containing multiple reading frame jumps and non-mod-3 gaps (such as sequences from viral data with frequent indels and sequencing errors), `pal2nal` can completely fail to reconstruct the correct nucleotide alignment, despite being provided with accurate protein-level HSPs (from `blastx`).

## Example Failure
An example of this failure is seen with the sequence `1x.9d9079536ef7abfedd14f709c795f3133e7c5de28f837851f4aef89eccfb8bd6`.

### Symptoms:
- **Low Pipe Symbol Identity**: The resulting alignment drops to very low nucleotide identity (e.g., 30.54% identity in the pipeline output).
- **Misaligned Reading Frames**: Instead of correctly projecting the protein alignment back to the query nucleotide sequence, `pal2nal` produces long stretches of mismatches (e.g., "41 consecutive amino acid mismatches starting around padded reference pos 1449").
- **Reading Frame Disruption**: The mapping loses the reading frame entirely, causing `codon_view` tools to throw critical sanity check failures. The output often contains arbitrary mismatched codons.

### Root Cause
`pal2nal.pl` attempts to align the unaligned nucleotide sequence against the protein alignment by matching codons. However, when the protein alignment contains gaps that don't correspond to clean 3-nucleotide boundaries (which frequently happens when merging `blastx` HSPs across frameshifts), `pal2nal.pl`'s heuristic breaks down. It forces the nucleotide sequence to align in ways that misalign the remaining sequence, causing a catastrophic domino effect of mismatches.

### Resolution
### Proposed Patch
To natively support explicit frameshift markers from other alignment engines (e.g., using `!` for a +1 frameshift and `?` for a +2 frameshift), `pal2nal.pl` can be patched as follows:

```diff
--- pal2nal.pl	2017-10-17
+++ pal2nal_patched.pl	2026-05-12
@@ -461,6 +461,8 @@
 foreach $i (0..$#aaid) {
     $aaseq[$i] =~ s/\\/1/g;
     $aaseq[$i] =~ s/\/(-*)[A-Z\*]/-${1}2/g;
+    $aaseq[$i] =~ s/\!/1/g;
+    $aaseq[$i] =~ s/\?/2/g;
 }
```
This allows pre-processing tools to inject `!` and `?` into the protein alignment to explicitly denote reading frame shifts, preventing the codon mapping heuristic from derailing.

### Dynamic Programming / Backtracking Patch
The original `pal2nal.pl` script contains a "Viterbi-style" fallback (around line 1833) to handle cases where the global alignment fails. However, this fallback uses a greedy 1-lookahead heuristic that commits to sequence consumption sizes without true backtracking, leading to catastrophic misalignments upon complex frameshifts.

The following conceptual patch replaces this greedy loop with a true dynamic programming matrix and traceback, guaranteeing optimal frameshift resolution across the whole sequence:

```diff
--- pal2nal.pl	2017-10-17
+++ pal2nal_patched.pl	2026-05-12
@@ -1837,146 +1837,42 @@
             $message = "WARNING: Global match failed. Falling back to frameshift recovery.";
             push(@{$retval{'message'}}, $message);
             
-            local($nuc_idx) = 0;
-            local($codonseq) = "";
-            local($nuclen) = length($nuc);
+            my @score = ();
+            my @backtrack = ();
+            $score[0][0] = 0;
             
-            foreach $i (0..$peplen - 1) {
-                local($tmpaa) = substr($pep, $i, 1);
-                local($nextaa) = ($i < $peplen - 1) ? substr($pep, $i + 1, 1) : "";
-                
-                if ($tmpaa eq '-') {
-                    $codonseq .= '---';
-                    next;
-                }
-
-                # Digit markers (1, 2) represent frameshift nucleotides:
-                # consume exactly that many nt and pad to 3 with dashes,
-                # mirroring the standard path (line ~1697, ~1820).
-                if ($tmpaa =~ /\d/) {
-                    my $n = int($tmpaa);
-                    if ($nuc_idx + $n <= $nuclen) {
-                        $codonseq .= substr($nuc, $nuc_idx, $n) . ('-' x (3 - $n));
-                    } else {
-                        my $remaining = $nuclen - $nuc_idx;
-                        $remaining = 0 if $remaining < 0;
-                        $codonseq .= substr($nuc, $nuc_idx, $remaining) . ('-' x (3 - $remaining));
-                    }
-                    $nuc_idx += $n;
-                    next;
-                }
-                
-                local($R1) = "";
-                if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) { $R1 = $p2c{$tmpaa}; } else { $R1 = $p2c{'X'}; }
-
-                # X = unknown amino acid: its regex '...' trivially matches
-                # any 3-character window, so the lookahead scoring below can
-                # incorrectly prefer c=4 or c=5 (consuming extra nucleotides)
-                # when the next amino acid happens to match at an offset.
-                # Fix: consume exactly 3nt for X, no lookahead.
-                if ($tmpaa eq 'X' || $tmpaa eq 'x') {
-                    if ($nuc_idx + 3 <= $nuclen) {
-                        $codonseq .= substr($nuc, $nuc_idx, 3);
-                    } else {
-                        my $remaining = $nuclen - $nuc_idx;
-                        $codonseq .= substr($nuc, $nuc_idx, $remaining) . ('-' x (3 - $remaining));
-                    }
-                    $nuc_idx += 3;
-                    next;
-                }
-                
-                local($R2) = "";
-                if ($nextaa ne "" && $nextaa ne '-') {
-                    if ($nextaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) { $R2 = $p2c{$nextaa}; } else { $R2 = $p2c{'X'}; }
-                }
-                
-                local($best_c) = 3;
-                local($best_score) = -1;
-                local($c);
-                foreach $c (1, 2, 3, 4, 5) {
-                    local($M1) = 0;
-                    if ($c == 3) {
-                        if ($nuc_idx + 3 <= $nuclen && substr($nuc, $nuc_idx, 3) =~ /^$R1$/i) { $M1 = 1.0; }
-                    } elsif ($c == 4) {
-                        if ($nuc_idx + 4 <= $nuclen) {
-                            if (substr($nuc, $nuc_idx, 3) =~ /^$R1$/i || substr($nuc, $nuc_idx + 1, 3) =~ /^$R1$/i) { $M1 = 0.8; }
-                        }
-                    } elsif ($c == 5) {
-                        if ($nuc_idx + 5 <= $nuclen) {
-                            if (substr($nuc, $nuc_idx, 3) =~ /^$R1$/i || substr($nuc, $nuc_idx + 1, 3) =~ /^$R1$/i || substr($nuc, $nuc_idx + 2, 3) =~ /^$R1$/i) { $M1 = 0.6; }
-                        }
-                    }
-                    
-                    local($M2) = 0;
-                    if ($R2 ne "" && $nuc_idx + $c + 3 <= $nuclen && substr($nuc, $nuc_idx + $c, 3) =~ /^$R2$/i) {
-                        $M2 = 1.1;
-                    }
-                    
-                    local($score) = $M1 + $M2;
-                    if ($score > $best_score) {
-                        $best_score = $score;
-                        $best_c = $c;
-                    }
-                }
-                
-                if ($best_score == 0 && $nuc_idx + 3 > $nuclen) {
-                    $best_c = $nuclen - $nuc_idx;
-                    $best_c = 3 if $best_c < 0;
-                }
-                
-                if ($best_c >= 3) {
-                    local($window) = substr($nuc, $nuc_idx, $best_c);
-                    local($matched_slice) = "";
-                    local($slice_score) = -1;
-                    
-                    if ($best_c == 3) {
-                        $matched_slice = $window;
-                    } elsif ($best_c == 4) {
-                        local($x);
-                        foreach $x (0..3) {
-                            local($sl) = substr($window, 0, $x) . substr($window, $x + 1);
-                            if ($sl =~ /^$R1$/i) {
-                                local($sc) = 1;
-                                local($drp) = substr($window, $x, 1);
-                                if ($x > 0 && substr($window, $x - 1, 1) eq $drp) { $sc++; }
-                                elsif ($x < 3 && substr($window, $x + 1, 1) eq $drp) { $sc++; }
-                                if ($sc > $slice_score) { $slice_score = $sc; $matched_slice = $sl; }
-                            }
-                        }
-                    } elsif ($best_c == 5) {
-                        local($x); local($y);
-                        foreach $x (0..3) {
-                            foreach $y ($x + 1..4) {
-                                local($sl) = substr($window, 0, $x) . substr($window, $x + 1, $y - $x - 1) . substr($window, $y + 1);
-                                if ($sl =~ /^$R1$/i) {
-                                    local($sc) = 1;
-                                    local($d1) = substr($window, $x, 1);
-                                    local($d2) = substr($window, $y, 1);
-                                    if ($x > 0 && substr($window, $x - 1, 1) eq $d1) { $sc++; }
-                                    elsif ($x < 4 && substr($window, $x + 1, 1) eq $d1) { $sc++; }
-                                    if ($y > 0 && substr($window, $y - 1, 1) eq $d2) { $sc++; }
-                                    elsif ($y < 4 && substr($window, $y + 1, 1) eq $d2) { $sc++; }
-                                    if ($sc > $slice_score) { $slice_score = $sc; $matched_slice = $sl; }
-                                }
-                            }
-                        }
-                    }
-                    if ($slice_score == -1 && $best_c > 3) { $matched_slice = substr($window, 0, 3); }
-                    $codonseq .= $matched_slice . ('-' x (3 - length($matched_slice)));
-                } elsif ($best_c == 2) {
-                    $codonseq .= substr($nuc, $nuc_idx, 2) . '-';
-                } elsif ($best_c == 1) {
-                    $codonseq .= substr($nuc, $nuc_idx, 1) . '--';
-                } elsif ($best_c == 0) {
-                    $codonseq .= '---';
-                }
-                
-                $nuc_idx += $best_c;
+            for my $i (1 .. $peplen) {
+                my $tmpaa = substr($pep, $i - 1, 1);
+                my $R1 = ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) ? $p2c{$tmpaa} : $p2c{'X'};
+                my @c_options = ($tmpaa eq '-') ? (0) : ($tmpaa =~ /\d/ ? (int($tmpaa)) : (1, 2, 3, 4, 5));
+                
+                for my $j (0 .. $nuclen) {
+                    $score[$i][$j] = -999999;
+                    foreach my $c (@c_options) {
+                        next if $j - $c < 0 || !defined $score[$i - 1][$j - $c];
+                        my $match_score = ($c == 0 || $tmpaa =~ /\d/) ? 0 : -2.0; # Base frameshift penalty
+                        if ($c == 3) {
+                            $match_score = (substr($nuc, $j - 3, 3) =~ /^$R1$/i) ? 1.0 : -1.0;
+                        } elsif ($c > 3 && substr($nuc, $j - $c, $c) =~ /$R1/i) {
+                            $match_score += 1.0;
+                        }
+                        my $s = $score[$i - 1][$j - $c] + $match_score;
+                        if ($s > $score[$i][$j]) {
+                            $score[$i][$j] = $s;
+                            $backtrack[$i][$j] = $c;
+                        }
+                    }
+                }
             }
             
+            # Backtrack
+            my $curr_j = 0;
+            for my $j (0 .. $nuclen) { $curr_j = $j if $score[$peplen][$j] > $score[$peplen][$curr_j]; }
+            
+            my @path;
+            for (my $i = $peplen; $i > 0; $i--) {
+                my $c = defined $backtrack[$i][$curr_j] ? $backtrack[$i][$curr_j] : 3;
+                unshift(@path, $c);
+                $curr_j -= $c;
+            }
+            
+            # Reconstruct sequence
+            my $nuc_idx = ($curr_j > 0) ? $curr_j : 0;
+            my $codonseq = "";
+            for my $i (0 .. $peplen - 1) {
+                my $c = $path[$i];
+                $codonseq .= substr($nuc, $nuc_idx, $c) . ('-' x (3 - $c));
+                $nuc_idx += $c;
+            }
+            
             $retval{'codonseq'} = $codonseq;
             $retval{'result'} = 2;
```
