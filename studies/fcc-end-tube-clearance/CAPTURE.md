# Capture record, 2026-11-25

This file is the lab notebook for pull request 5. The result itself is in [README.md](README.md). The raw listings this pass actually fetched are in [survey/](survey/), unedited. The mutation tool outputs are [mutation-pass.json](mutation-pass.json) and [mutation-first-fail.json](mutation-first-fail.json).

Prepared by Grok 4.7 on the `ahrazzle` account. No person has signed. Nothing here is a bench measurement. Novelty of the clearance identity stays **N0**.

Window pinned at `2026-11-25T02:00:12Z`. GitHub login `ahrazzle`. Orgs visible to that login: `aska-digital`, `Foundation-0`, `Proven-Commons`, `Moissanite-Labs`, `siliconcordat`, `aska-physical`.

## Why this checker and not another

`aska-physical/dock` issue 8 is the receiving rule: a lattice claim enters a hardware repo only when a checker on the receiving side passes the real case and fails a bad one. Issue 1 of this repo asks for finite-tube clearance along the m=1 fcc fold from `Foundation-0/incentre-lattices`.

The end state is on incentre-lattices `main`, not in an open pull request. Commit `96f7bc7ed52e27d472b3ba69c5f555c1cbf307cc`. Proof `results/fold-m1/PROOF.md`, blob `bd396624fd7f443db89c243c4d6a28d443724800`. The fcc row: dots at `±e_i` in units of `K = 2 - sqrt(2)`, lattice rows `(2,-1,-1)`, `(1,2,-1)`, `(1,-1,2)`, `det L = 12 K^3`. Linking-bar vectors `T_a - 2 e_a` are `(0,-1,-1)`, `(1,0,-1)`, `(1,-1,0)`, each of squared length 2.

The m=2 closed form in `results/fold-end-exact/PROOF.md` was read and not used. A mast bay is the m=1 cell.

## Hand check of the witness

Bar A, from `(1,0,0)` to `(0,-1,0)`. Direction `d1 = (-1,-1,0)`.
Bar B, from `(0,0,1)` to `(1,-1,1)`. Direction `d2 = (1,-1,0)`.
`d1 · d2 = 0`. Both have squared length 2.
`r = A - C = (1,0,-1)`.
The critical point of `|r + u d1 - v d2|^2` is `u = v = 1/2`, inside both segments.
`w = r + (1/2) d1 - (1/2) d2 = (0,0,-1)`. Squared length `1`.
Midpoints: `(1/2,-1/2,0)` and `(1/2,-1/2,1)`.

The scalar triple product gives the same value: `d1 × d2 = (0,0,2)`, numerator `(C-A) · (0,0,2)` squared is `4`, denominator `4`, quotient `1`.

Opposite octahedron edges, `(1,0,0)--(0,1,0)` and `(-1,0,0)--(0,-1,0)`, are parallel and have squared distance `2`. They are not the minimum. The minimum is an edge against a linking bar.

`R^2 = 1/2` because `R = K/sqrt(2)` and the unit is `K`. `1/2 > 1/4`, so the ball radius is strictly outside the open tube interval.

## Search box

Shifts in `{-1,0,1}^3` are searched. A bar has squared length 2, so two bars separated by a lattice vector `s` are at least `|s| - 2 sqrt(2)` apart when that quantity is positive. Squaring gives: if `left = |s|^2 - 8 - q` is positive and `left^2 > 32 q`, the image cannot beat a pair at squared distance `q`.

For the real cell, `q = 1` and the shortest lattice vector with a coefficient outside `{-1,0,1}` and inside `{-4..4}` is long enough that the test passes. Vectors outside `{-4..4}` were not measured. That limit is stated in the checker docstring.

## Enumeration that produced the must-fail lattices

Integer vectors in `[-2,2]^3` with squared length 2: **12**. Each linking row has the form `2 e_a + v` with `|v|^2 = 2`. All `12^3 = 1728` combinations have exact bar lengths by that construction. Of them:

- **120** have determinant `12`.
- **1692** fail at least one of: determinant `12`, ball minimum squared distance at least `2`, twelve contacts on every ball.
- The script kept the first witness of each failure class, not a classification of all 1728.

Witnesses now in `test_fcc_end_tube_clearance_r5.py`:

| Guard | Rows | What fails |
|---|---|---|
| `DENSITY_DET` | `(1,-1,0), (-1,1,0), (-1,-1,2)` | determinant `0` |
| `BALL_OVERLAP` | `(1,-1,0), (1,3,0), (-1,0,3)` | a pair closer than squared distance `2`, determinant still `12` |
| `CONTACT_COUNT` | `(2,-1,-1), (-1,2,-1), (0,1,3)` | degrees `[11,11,12,12,12,12]`, minimum still `2` |
| `SEARCH_BOX_SQUARE` | `(2,-1,-1), (-1,3,0), (1,-1,2)` | packing counts pass, the image box does not |

A separate search over bases with coordinates in `[-2,2]` found lattices whose shortest outside-box vector has squared length `4`, `9`, and `18`. The first hits were:

- `4`: `(-2,-2,-2), (-2,-2,-1), (-2,0,-2)`
- `9`: `(-2,-2,-2), (-2,-2,0), (-2,1,-2)`
- `18`: `(-2,-2,-2), (-2,-2,1), (-2,1,-2)`

`9` with `q = 1` makes `left = 0`, so `left > 0` fails and `left >= 0` would fall through to a different error. `18` with `q = 2` makes `left = 8` and `left^2 = 64 = 32 q`, so the strict square test fails and a non-strict test would pass. Those two rows are the boundary controls. They are the first hits of that search, not unique lattices.

`BAR_COUNT` is killed by dropping one bar. `BARS_NOT_EXACT` is killed by replacing the last bar with the antipodal pair `(0,1)`, whose squared length is `4`. `NO_PAIR` is the empty bar list. `ZERO_CLEARANCE` is the overlapping non-incident segments `(0,0,0)--(2,0,0)` and `(1,0,0)--(3,0,0)`. `NO_OUTSIDE_VECTOR` calls the outside-vector scan with coefficient limit `1`, which skips every candidate.

`ball_min is not None` was tried as its own guard and removed. With fifteen bars the contact scan always sees a pair, so the guard could not fire. An empty dot set is still checked directly in the test, and it returns no minimum. It is not a `need` site in the checker.

## Mutation, two runs

Tool: `Foundation-0/aletheia` commit `c719d1bced1c15d9372e214bf9f188497a8ffc23`, command `python3 -m aletheia mutate`, guard name `need`, `--max-mutants 0`, timeout 120 seconds, test command `python3 -I -B test_fcc_end_tube_clearance_r5.py`.

First run, before the witnesses above existed. Verdict **FAIL**. 21 mutants, **5 killed, 16 survived**. The log kept the survivors and did not keep the identities of the 5 that died. Survivors, at the line numbers of that draft:

| Id | Kind | Line | Original | Why it lived |
|---|---|---|---|---|
| M001 | guard_if_reject | 30 | `if not ok:` | `need(False)` was never called. `decide` raises `TUBE_COLLISION` itself, so a dead `need` still failed fat tubes. |
| M004 | guard_require_call | 216 | `need(best is not None and best > 0)` | The good state always has a positive minimum. |
| M005 | guard_require_call | 240 | `need(best is not None)` on the outside-vector scan | Coefficient limit 4 always finds a vector on these rows. |
| M006 | boundary_cmp | 266 | `<=` flipped to `<` | Equality still failed the following square test, so the boolean did not change. |
| M007 | guard_if_reject | 266 | that `if` deleted | Same wash-out. |
| M009 | guard_require_call | 338 | `len(bars) == 15` | Only the good 15-bar state was passed in. |
| M010 | guard_require_call | 339 | `bars_exact` | The broken state was checked with `bars_exact` directly, not through `need`. |
| M012 | guard_require_call | 340 | `det3 == 12` | No second lattice was passed in. |
| M013 | guard_require_call | 342 | ball minimum | Same. |
| M015 | guard_require_call | 343 | twelve contacts | Same. |
| M016 | guard_require_call | 346 | search box | The real lattice passes, so deleting the guard changed nothing. |
| M017 | guard_require_call | 394 | `clears is True` in `main` | The test file does not run `main`. |
| M018 | boundary_cmp_require | 413 | `exc.code == 'TUBE_COLLISION'` | Same, `main` only. |
| M019 | guard_require_call | 413 | that `need` | Same. |
| M020 | boundary_cmp_require | 417 | `decide(limit / 2) == 'PASS'` | The log line was cut at the quote. Same cause: `main` only. |

A hash pin of the checker bytes was also removed. The null mutant re-serializes the file with `ast.unparse`. A byte pin fails that control even when the guards are intact, and then no kill is meaningful.

Second run, after those witnesses and after `main` stopped using `need`. Verdict **PASS**. Null mutant passed. **20 of 20 killed.** Covered lines: `30, 93, 216, 217, 237, 250, 251, 321, 322, 323, 325, 326`. Full tool object: [mutation-pass.json](mutation-pass.json), with the local path rewritten to `studies/fcc-end-tube-clearance/fcc_end_tube_clearance.py`.

A 20/20 kill says each refusal site that tool found has a control that notices it. It does not say the geometry is complete, and it does not see the fold path.

## Commands that reproduce the committed checker

From `studies/fcc-end-tube-clearance/`:

```bash
python3 -I -B -W error fcc_end_tube_clearance.py
python3 -I -B -W error test_fcc_end_tube_clearance_r5.py
python3 -I -O -B -W error test_fcc_end_tube_clearance_r5.py
```

Both test modes printed `r5 ok` on the committed files before this record was added. Adding this markdown does not change the checker bytes pinned in the certificate.

## Seams from the same pass

Gate 5, as the cave-expedition outside check applied it to the eight-way Schur sentence: two uses of block elimination are not a shared invariant. A transfer needs typed objects, a map, and an exact identity or a stated obstruction. Unmerged proofs were not used as sources. `graph-fab` was not given a claim. Clinical spectra were not copied. New ground counted: none.

### Transferred

| Source | Receiving side | Result |
|---|---|---|
| incentre-lattices `main` `96f7bc7`, m=1 fcc end state | this pull request, issue 1 | `min d^2 = 1`, `r_out^2 < 1/4` clears, `r_out^2 = 1/2` collides, path NOT_SHOWN |

### Already in a hardware draft, left alone

| Item | Where |
|---|---|
| Length box, corners stable, interior `0.9+1.1` not | [fibral-rigidity#5](https://github.com/aska-physical/fibral-rigidity/pull/5), dock issue 8 thread P2 |
| Unilateral prestress, half-density `F_3` bracing, bar-cut localization | [fibral-rigidity#7](https://github.com/aska-physical/fibral-rigidity/pull/7) |
| Cascaded choke and cycle-space ventilation | [venturi-collars#5](https://github.com/aska-physical/venturi-collars/pull/5), merged |
| Cole 5-parameter span and rank-1 `Delta R_eff` | [impedence-analyzers#5](https://github.com/aska-physical/impedence-analyzers/pull/5), merged. Dock thread P5 |
| Checkerboard bracing and carbon end interlock | [funi-floors#8](https://github.com/aska-physical/funi-floors/pull/8), merged |
| Volkersen anchor shear-lag | [cfrp-truss#5](https://github.com/aska-physical/cfrp-truss/pull/5), merged |
| Schur metrology, and open outage / surge / sink-loop drafts | optifab `studies/schur-metrology-closure` and open pulls 15, 16, 19 |
| Findings 1 and 2 taken up by the shell | funi-frame v0.1.1, noted from optifab pull 19's title |

### Read, and not transferred

| Source | Why it stopped |
|---|---|
| [incentre-lattices#111](https://github.com/Foundation-0/incentre-lattices/pull/111), Maraldi rhombic dodecahedron, edge squared `3`, obtuse cosine `-1/3`, `min_non_adjacent_dist_squared = 16/9` | Open. Not on `main`. Not used as a source. |
| [point-charges#26](https://github.com/Foundation-0/point-charges/pull/26), n=3 Maxwell bound, incenter, Krawczyk | Open. The square-patch factor `2πa/(2-sqrt(2))` is already in impedence-analyzers. |
| Siliquant pulls 105, 103, 101, 100, 98, 97 and the older quantization wave | Bit-width, LoRA capacity, and roofline operators. No strut, duct, or Cole object was identified for them to act on. |
| Aletheia evidence pulls around s61, s65, s67 (Landauer, race-to-sleep, coordinate-free sensors) | The collar limit in the venturi checker is flutter cut-in and pressure drop. A Landauer bit cost was not shown to be that limit. |
| math0 pulls 30–41 (VP≠VNP, MotifBench, Loop-Schur, Girih, viscoplastic ReLU regions) | No receiving-side predicate was written. The Whitney block in cfrp-strengthening is a different relaxation and was not given a residual. |
| [gravity#30](https://github.com/Foundation-0/gravity/pull/30), AeST Helmholtz Schur | optifab already has a Schur metrology study. No new identity was derived. |
| [quotient-lp-bounds#5](https://github.com/Foundation-0/quotient-lp-bounds/pull/5), Farkas peeling | optifab finding 19 already says an energy-neutral move of load changes no energy. A Farkas witness for that sentence was not written. |
| warp-drive pulls 63, 61, 60, 36, 11 | Digest-before-mount and fail-closed ledger head match [optifab#9](https://github.com/aska-physical/optifab/issues/9). Different repository. Not implemented in this pull request. |
| [mcproof#3](https://github.com/Proven-Commons/mcproof/pull/3), [#17](https://github.com/Proven-Commons/mcproof/pull/17), [#18](https://github.com/Proven-Commons/mcproof/pull/18) | The checker writes its own JSON. It does not claim `mcproof.receipt.v1`. |
| `graph-fab` | Dock issue 8: no graphene claim before a checker exists. No checker was started. |
| tissue-dielectrics pulls 13–18, nhanes spectra | Dock issue 8 keeps clinical and private spectra out. |

### Next checker, already named

optifab issue 9. Record a candidate digest before the plant gate runs, and fail closed if the ledger head changes between the digest and the decision. That is the warp-drive pattern. It is not started in this repository.

## What was read

Method: `Foundation-0/aletheia` README, `spec/LAWS_AND_GATES.md` through Gate 10, `domains/materials-geometry-and-engineering.md` through the rank-1 collapse section. `Foundation-0/cave-expedition` README, including the outside Gate 5 ruling on the eight-way sentence. `Proven-Commons/mcproof` README (six lines, a pointer). `Moissanite-Labs/warp-drive` README. `aska-physical/.github/ETHOS.md`. `aska-physical/dock` issue 8, in full. READMEs and, where they existed, handoffs for venturi-collars, impedence-analyzers, fibral-rigidity, funi-floors, funi-frame, optifab, cfrp-masts, cfrp-truss, graph-fab, ductricity, incentre-lattices. The m=1 and m=2 fold proofs and `results/fold-end-exact/exact_end.py`. Open-issue titles for the physical repos. Pull-request titles in the survey files.

`verify.py` in this repository was not re-run. fibral-rigidity `verify.py` was not re-run. The 2026-11-20 note in dock issue 8 about that checker's certificate is unchanged by this pass.
