# m=1 fcc end state: non-incident tube clearance

Receiving-side checker for [cfrp-masts#1](https://github.com/aska-physical/cfrp-masts/issues/1), under the rule in [aska-physical/dock#8](https://github.com/aska-physical/dock/issues/8). A lattice claim crosses into this repository when a checker here passes the real case and fails a bad one. This is a model. It is not a bench measurement. Nothing in it is a product radius in millimetres.

Novelty label: **N0**. Oracle: **exact**. The closed form is the fcc row of `Foundation-0/incentre-lattices` `results/fold-m1/PROOF.md` at commit `96f7bc7ed52e27d472b3ba69c5f555c1cbf307cc` (blob `bd396624fd7f443db89c243c4d6a28d443724800`). Segment distance is classical. The local fact is a scope correction: the proved ball packing does not clear a tube whose outer radius equals the proved ball radius.

## The state

One octahedron. Units of `K = 2 - sqrt(2)`. Dots at `±e_i`. Lattice rows `(2,-1,-1)`, `(1,2,-1)`, `(1,-1,2)`. `det L = 12 K^3`. Fifteen bars: the twelve octahedron edges, and three linking bars from `+e_a` to `-e_a` in the next cell. Every bar has squared length `2`, which is `(2R)^2` for `R = sqrt(2) - 1`. The checker recomputes that, the determinant, and twelve contacts per ball with no closer pair. It does not read the upstream repository.

## What is certified

Non-incident bars do not share a vertex. Shared vertices are joints; their centerline distance is `0`, and counting them would make every tube look like a collision.

Over the fifteen bars and lattice shifts in `{-1,0,1}^3`, the minimum squared centerline distance is **`1`**. Six unordered pairs attain it. One witness, checked by hand:

- edge from `(1,0,0)` to `(0,-1,0)`
- linking bar from `(0,0,1)` to `(1,-1,1)`
- midpoints `(1/2, -1/2, 0)` and `(1/2, -1/2, 1)`
- separation `(0,0,-1)`

So the centerlines come **one `K` apart**. They meet when each outer radius reaches `K/2`, that is `r_out^2 = 1/4`. The proved balls have `R^2 = 1/2` in the same units, and `1/2 > 1/4`, so a tube as thick as those balls meets the other strut. A tube with `r_out^2 < 1/4` does not, at this end state.

`r_out / L_bar < sqrt(2)/4`, strict, because `L_bar^2 = 2` and `(2 r_out)^2 < 1`.

The search box is exhaustive for this minimum: every lattice vector with a coefficient outside `{-1,0,1}` and inside `{-4..4}` is long enough that a bar carried by it stays farther than `1` from a fundamental bar, by the triangle bound `|s| - 2 sqrt(2)`.

## Pass and fail

| Case | Result |
|---|---|
| `r_out^2 = 0` | PASS |
| `r_out^2 = 1/5` | PASS |
| `r_out^2 = 1/4` (the witness) | FAIL `TUBE_COLLISION` |
| `r_out^2 = 1/2` (the proved ball) | FAIL `TUBE_COLLISION` |
| one dot moved by `1` | FAIL `BARS_NOT_EXACT` |

A second distance formula (scalar triple product, written in the test file) agrees on every non-incident pair. Opposite octahedron edges stay at squared distance `2`. They are not the witness. The witness is an edge against a linking bar, so the tight gap is mid-span, not at the hinge.

## What this does not show

- The path. `fold-m1` certifies ball centers along a sampled trajectory. It does not certify segment distances. The interior of the fold is **NOT_SHOWN**.
- A collar at a shared vertex. That clearance is identically zero; it is a joint, and this checker leaves it out.
- Any radius in the product checker. `verify.py` is not edited and was not re-run. No millimetre figure is derived here.
- New ground. Dock issue 8 counts new ground only after a merge and a second reviewer. This label stays `N0`.

## Examined and not transferred

The same pass looked at the live foundation-0, moissanite-labs, and proven-commons streams (open and merged, through 2026-11-25). These were not given a typed map into this checker. Recording them here is so the next agent does not re-derive the refusal.

- `incentre-lattices` pull 111 (Maraldi rhombic dodecahedron) is open. An unmerged proof was not used as a source.
- `fibral-rigidity` pulls 5 and 7 already carry the length-box and the prestress certificate. Not copied.
- `point-charges` pull 26 (Maxwell bound, Krawczyk) is open. Not used.
- Siliquant quantization studies, and the aletheia silicon Landauer and race-to-sleep tracks, have no strut or bar object in this state. No operator was written. Gate 5 of the cave-expedition protocol refuses a shared word as a bridge. The outside check on the eight-way Schur sentence is the precedent: block elimination in two places is not a common invariant.
- `warp-drive` fail-closed digest-before-mount and `mcproof` receipts match `optifab` issue 9, which is a different repository. This file does not claim `mcproof.receipt.v1`.
- `graph-fab` stays untouched. Dock issue 8 says no graphene claim before a checker exists.

## Mutation

`Foundation-0/aletheia` at `c719d1bced1c15d9372e214bf9f188497a8ffc23`, command `aletheia mutate` with guard name `need`, `--max-mutants 0`. The null mutant passed. 20 of 20 mutants were killed. That says each refusal site the tool found has a control that notices it. It does not say the geometry is complete.

## Reproduce

From this directory, Python 3, standard library only:

```bash
python3 -I -B -W error fcc_end_tube_clearance.py
python3 -I -B -W error test_fcc_end_tube_clearance_r5.py
python3 -I -O -B -W error test_fcc_end_tube_clearance_r5.py
```

`need` raises `CheckFailure`. `python3 -O` does not remove it.

## Credit

Prepared by Grok 4.7, acting on the `ahrazzle` account, for the aska-physical receiving rule. No person has signed this. Human sign-off on the mast repository remains Shahid Arifuddin (`@abu-aska`). The closed form belongs to `Foundation-0/incentre-lattices`.
