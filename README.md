# cfrp-masts

**Snap-Locked Deployable CFRP Masts, Telecom Towers & Space Booms.**

Translates the exact continuous folding geometry of [`Foundation-0/incentre-lattices`](https://github.com/Foundation-0/incentre-lattices), the area-graded tubular buckling & octagonal-pore cable laws of [`aska-physical/fibral-rigidity`](https://github.com/aska-physical/fibral-rigidity), and the rational linkage state certification of [`Foundation-0/rigid-realizations`](https://github.com/Foundation-0/rigid-realizations) into full-scale deployable structural booms and masts.

## Why Existing Market Products Fall Short

Current lightweight mast and space-boom technologies force a sharp trade-off:
1. **Fixed filament-wound composite trusses** ([Coreless Filament-Wound Space Trusses](https://pmc.ncbi.nlm.nih.gov/articles/PMC11001885/)) and **3D-printed metal lattices** ([Deshpande, Fleck & Ashby 2001](https://people.bu.edu/parkhs/Papers/heJMPS2017.pdf)) cannot fold compactly (`1.0x` compaction) and carry heavy metallic mass or thermal expansion drift.
2. **Slit-tube (STEM), collapsible-tube, and inflatable booms** ([Deployable Boom Benchmarks](https://matheo.uliege.be/bitstream/2268.2/5447/5/Thesis_Zain_Mansoor.pdf)) stow small but suffer low torsional/shear stiffness, deployment hysteresis, and millimeter-to-centimeter thermal pointing drift.

## Verified Product Architecture & Market Superiority

Each mast bay uses thin-walled pultruded/filament-wound CFRP tubes (inner/outer radius ratio $\eta = 0.80 \implies 4.56\times$ higher local Euler buckling resistance than solid rods) graded at $\alpha_K^* = 2\sqrt{2}$ with 12 non-intersecting $45^\circ$ octagonal-pore [Torayca T1100G](https://www.toraycma.com/wp-content/uploads/T1100G-Data-Sheet.pdf) or [Mitsubishi DIALEAD K13D2U](https://www.m-chemical.co.jp/carbon-fiber/en/product/tow/) cables ($\beta_{\text{iso}} = 0.362323$) that fold slack during stowage and **snap taut at full deployment** to lock shear rigidity ($GA_{\text{eff}} > 0$).

### Verified `5.0 m` Space Magnetometer / Solar Boom Benchmark ($0.25\text{ m}$ Cross-Section, $\bar{\rho} = 2.5\%$, $5\text{ kg}$ Tip Payload, $50\text{ N}$ Slew Load)

| Metric | Al 6061-T6 Octet Boom | 3D-Printed Ti-6Al-4V Octet | `cfrp-masts` K13D2U Snap-Lock Boom | Verified Advantage |
| :--- | :---: | :---: | :---: | :--- |
| **Structural Mass** | $21.09\text{ kg}$ | $34.61\text{ kg}$ | **$14.06\text{ kg}$** | **59% lighter than Ti** |
| **Flexural Rigidity $EI$** | $62.3\text{ kN}\cdot\text{m}^2$ | $102.9\text{ kN}\cdot\text{m}^2$ | **$586.1\text{ kN}\cdot\text{m}^2$** | **14.02x specific $EI/\text{kg}$ vs. Ti** |
| **Tip Deflection ($50\text{ N}$)** | $33.47\text{ mm}$ | $20.27\text{ mm}$ | **$3.56\text{ mm}$** (Timoshenko exact) | **5.7x–9.4x stiffer** |
| **1st Natural Frequency $f_1$** | $1.94\text{ Hz}$ | $2.17\text{ Hz}$ | **$6.52\text{ Hz}$** | **3.0x higher control bandwidth** |
| **Stowed Stack Length** | $5.00\text{ m}$ (fixed) | $5.00\text{ m}$ (fixed) | **$1.51\text{ m}$** ($3.32\times$ fold) | **Fits compact fairing / truck bed** |
| **$100\text{ K}$ Thermal Drift** | $11.55\text{ mm}$ | $4.30\text{ mm}$ | **$0.30\text{ mm}$** | **38.5x lower pointing drift vs. Al** |

`verify.py` also sizes and verifies a **`15.0 m` portable T1100G telecom/defense mast** ($0.60\text{ m}$ bay, $251.7\text{ kg}$, stowing to **$4.52\text{ m}$**, global Haringx buckling capacity **$88.99\text{ kN}$**).

## Reproduce & Verify

```bash
python3 verify.py
```

## `aska-physical` Ecosystem & Upstream Lineage

| Repository | Role in Physical Stack | Upstream `Foundation-0` Lineage |
| :--- | :--- | :--- |
| [`aska-physical/fibral-rigidity`](https://github.com/aska-physical/fibral-rigidity) | Continuum CFRP truss homogenization, area-grading (`alpha_K* = 2*sqrt(2)`, `alpha_E* = 2*sqrt(6)/3`), and 45-deg pore-tie theory | [`Foundation-0/incentre-lattices`](https://github.com/Foundation-0/incentre-lattices) |
| [`aska-physical/cfrp-strengthening`](https://github.com/aska-physical/cfrp-strengthening) | Step-by-step re-runnable ACI 440.2R concrete beam strengthening baseline | [`Foundation-0/cfrp-strengthening`](https://github.com/Foundation-0/cfrp-strengthening) |
| [`aska-physical/cfrp-truss`](https://github.com/aska-physical/cfrp-truss) | Truss-anchored NSM + 45-deg shear-stitched, self-sensing CFRP beam/column retrofit kits | `cfrp-strengthening` + `fibral-rigidity` + `nhanes-bioimpedance` |
| [`aska-physical/cfrp-masts`](https://github.com/aska-physical/cfrp-masts) | Snap-locked 3.32x-compacting deployable CFRP booms, telecom towers, and solar masts | `incentre-lattices` + `fibral-rigidity` + [`Foundation-0/rigid-realizations`](https://github.com/Foundation-0/rigid-realizations) |
| [`aska-physical/venturi-collars`](https://github.com/aska-physical/venturi-collars) | Battery-free aeroelastic M55J/PCB-coil Venturi collars for VAV HVAC & industrial ducts | [`Foundation-0/ductricity`](https://github.com/Foundation-0/ductricity) + [`Foundation-0/thermal-capture`](https://github.com/Foundation-0/thermal-capture) |
| [`aska-physical/impedence-analyzers`](https://github.com/aska-physical/impedence-analyzers) | Artifact-rejecting 7-frequency Cole EIS analyzers for structural NDT & clinical BIA | [`Foundation-0/nhanes-bioimpedance`](https://github.com/Foundation-0/nhanes-bioimpedance) + [`Foundation-0/point-charges`](https://github.com/Foundation-0/point-charges) |

## Attribution & Ethos

- **Architect & Human Sign-Off:** Shahid Arifuddin (`@abu-aska`)
- **Foundational `f0` Lineage:** Ahraz Arifuddin (`@ahrazzle`)
- **Autonomous Agent Capture:** Google Gemini (per Aska inter-agent collaboration ethos; see [`HANDOFF.md`](HANDOFF.md) and [`CHANGELOG.md`](CHANGELOG.md)).
