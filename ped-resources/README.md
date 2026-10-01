# PED resources

Ready-to-run notebooks for analysing ensembles from the [Protein Ensemble Database (PED)](https://proteinensemble.org) with SOURSOP. They run in Google Colab with nothing installed: click a badge, set the PED identifier, and run.

| Notebook | Colab | What it does |
| --- | --- | --- |
| `ped_ensemble_characterisation.ipynb` | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/ped-resources/ped_ensemble_characterisation.ipynb) | A standard biophysical characterisation of any PED ensemble, chosen by its identifier, aimed at users who have never written code (all code is hidden). A metapredict disorder prediction (with a warning if any region is predicted folded); Rg and end-to-end distance distributions against length-matched SAW-ν models (ν = 0.33, 0.5, 0.589); per-residue helicity and β-strand content from DSSP; the distance map and a scaling map normalised by the AFRC, each with a slider for the colour scale; the distance between any two residues (chosen with sliders) against the AFRC; and Rg against asphericity over a large AFRC ensemble. Uses SOURSOP, [afrc](https://github.com/idptools/afrc) and [metapredict](https://github.com/idptools/metapredict). |

For a tour of the PED API itself (searching, loading, weights and so on), see the [PED examples](../demo_examples/ped_examples/).
