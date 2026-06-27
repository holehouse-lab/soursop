# SAXS reweighting of a STARLING ensemble (PNt)

Worked examples of deriving per-conformer **weights** for a disordered-protein ensemble
from small-angle X-ray scattering (SAXS) data, using SOURSOP's `ssbme` reweighters. The
system is the N-terminal region of PNt; the conformational ensemble was generated with
STARLING.

Of note, [STARLING](https://idptools-starling.readthedocs.io/en/latest/) enable the generation of coarse-grained conformational ensembles of IDRs based on just sequence as input. While it's generally quite accurate, when experimental agreement is sub-par we can use BME to reweight ensembles to further improve agreement (see Figure 5 in the STARLING paper) .

STARLING provides BME built in as direct support, but does not provide iBME, BMECustom or
COPER as implemented here in SOURSOP.

These notebooks reweight the *same* ensemble against the *same* SAXS measurement in several
different ways, so you can see what changes when you fit the full scattering curve versus a
single derived radius of gyration, and how the different reweighting algorithms compare.

## Notebooks

| notebook | reweighter | fits against | take-away |
|----------|------------|--------------|-----------|
| [`bmecustom_saxs_reweighting.ipynb`](bmecustom_saxs_reweighting.ipynb) | `BMECustom` | the **full raw SAXS curve** (per-frame profiles) | how to plug a custom, scale-marginalised goodness-of-fit into `BMECustom`; weights then flow into any SOURSOP observable |
| [`bme_rg_reweighting.ipynb`](bme_rg_reweighting.ipynb) | `BME` | a **single derived Rg** (Guinier or MFF) | how standard `BME` reweights to a scalar observable, and how that differs from using the raw curve |
| [`coper_saxs_reweighting.ipynb`](coper_saxs_reweighting.ipynb) | `COPER` | the **full raw SAXS curve** (pre-scaled) | the hard-`χ²`-constraint alternative to BME: pick a target fit, get the highest-entropy ensemble that reaches it, and a feasibility test |
| [`reweighting_method_comparison.ipynb`](reweighting_method_comparison.ipynb) | `BME`, `iBME`, `BMECustom`, `COPER` | the **full raw SAXS curve** | the four reweighters head-to-head at matched effective sample size: do they produce the same ensemble? |

Run them top-to-bottom; each is self-contained and reads the data files in this directory.

### Headline comparison 1 — single Rg vs the raw curve (from `bme_rg_reweighting.ipynb`)

| weighting | ⟨Rg⟩ (Å) | phi | full-curve χ² |
|-----------|---------:|----:|--------------:|
| prior (uniform) | 57.1 | 1.00 | 4.40 |
| BME → Guinier Rg | 49.9 | 0.81 | 2.94 |
| BME → MFF Rg | 51.1 | 0.87 | 3.05 |
| BMECustom (full curve) | 52.9 | 0.77 | **1.13** |

Matching a single Rg pulls the mean size onto the target but only partly improves agreement
with the measured scattering; fitting the raw curve uses strictly more information, reaches a
much better χ², and settles at a different, data-supported ⟨Rg⟩.

### Headline comparison 2 — the four reweighters at matched φ ≈ 0.6 (from `reweighting_method_comparison.ipynb`)

| method | scale handling | ⟨Rg⟩ (Å) | full-curve χ² |
|--------|----------------|---------:|--------------:|
| BME (prescaled) | fixed | 56.2 | 1.39 |
| COPER (prescaled) | fixed | 56.2 | 1.39 |
| iBME | fitted | 58.3 | 1.57 |
| BMECustom | marginalised | 51.4 | **0.76** |

At a matched operating point **BME and COPER give the same ensemble** (weight correlation
≈ 1.0) — the penalty and hard-constraint forms trace one maximum-entropy curve. The visible
differences come from **how the SAXS scale is treated** (fixed vs fitted vs marginalised),
not from the choice of reweighting algorithm.

## Data files

| file / folder | description |
|---------------|-------------|
| `pnt.dat` | **Experimental** SAXS profile — three columns: `q` (Å⁻¹), `I(q)`, `σ(I)`. |
| `pnt_STARLING.pdb`, `pnt_STARLING.xtc` | The STARLING conformational ensemble (600 conformers, topology + trajectory). |
| `saxs_frames/frame_NNNNN.dat` | One **computed** scattering curve per conformer (600 files), two columns: `q`, `I(q)`. `frame_00000.dat` ↔ trajectory frame 0, etc. Stacked into the `(600, m)` matrix `BMECustom` expects. |
| `pnt_starling_average_scattering.dat` | The uniform-weight ensemble-average computed profile (`q`, `I(q)`, `σ`), for reference. |
| `guiner_analysis/pnt_guinier.txt` | Guinier-analysis report; the notebooks parse `Rg = 49.85 ± 0.37 Å` from it. (`.png`/`.pdf` are the analysis plots.) |
| `mff_analysis/pnt_MFF_fit.txt` | Molecular-form-factor fit report; the notebooks parse `Rg = 51.11 ± 0.13 Å` from it. |

> The computed and experimental curves are on different `q`-grids and have an arbitrary
> overall intensity scale, so the notebooks interpolate onto a common grid and marginalise the
> scale inside the cost function.

## Requirements

`soursop`, plus `numpy` and `matplotlib`. The notebooks were executed with these installed; to
re-run, open them in Jupyter (or `jupyter nbconvert --to notebook --execute`).

## Key SOURSOP entry points used

- `soursop.ssbme.BMECustom` — reweight against a vector/matrix observable with a pluggable cost.
- `soursop.ssbme.BME` / `soursop.ssbme.iBME` — penalty-form reweighting against scalar
  observables; `iBME` additionally fits the SAXS scale/offset.
- `soursop.sscoper.COPER` / `soursop.sscoper.iCOPER` — hard-`χ²`-constraint reweighting (with a
  feasibility test); `iCOPER` fits the scale/offset.
- `soursop.ssbme.ExperimentalObservable` — the shared (value, uncertainty) observable type.
- `soursop.sstrajectory.SSTrajectory` — load the ensemble; `get_radius_of_gyration()` returns
  the per-frame Rg, and accepts `weights=` to return the reweighted ensemble average.

See the SOURSOP docs (`modules/bme`, `modules/coper`) for the full reweighting API.
