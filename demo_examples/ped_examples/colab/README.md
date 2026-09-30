# PED examples for Google Colab

#### Last updated 2026-09-30

These PED example notebooks are set up to run in Google Colab, so you can try them without installing anything. Click a badge to open a notebook in Colab, then run it from the top. The first code cell installs SOURSOP with [uv](https://docs.astral.sh/uv/), which usually takes well under a minute.

| Notebook | |
| --- | --- |
| 1. Getting started | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/01_getting_started.ipynb) |
| 2. Searching PED | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/02_searching_ped.ipynb) |
| 3. One protein, many methods | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/03_one_protein_many_methods.ipynb) |
| 4. Chain dimensions across PED | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/04_chain_dimensions_across_ped.ipynb) |
| 5. Distance and contact maps | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/05_distance_and_contact_maps.ipynb) |
| 6. Local structure and NMR observables | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/06_local_structure_and_nmr.ipynb) |
| 7. Weighted ensembles and working with files | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/holehouse-lab/soursop/blob/master/demo_examples/ped_examples/colab/07_weighted_ensembles_and_downloads.ipynb) |

Each Colab session starts fresh, so the notebooks download their ensembles again every time. Most take under a minute; notebook 4 fetches 30 ensembles and takes around 5 minutes. The notebooks are the same as the local versions, except for the install cell. To see what the output should look like, open the local versions one folder up, which are saved with their output.

These notebooks are generated from the local ones by `make_colab_notebooks.py`, so edit the local notebooks and then run `python make_colab_notebooks.py` here, rather than editing these directly.
