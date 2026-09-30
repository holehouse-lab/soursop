ped - the Protein Ensemble Database
=========================================================

Overview
----------------------------

The `Protein Ensemble Database (PED) <https://proteinensemble.org>`_ is the main public repository of conformational ensembles of intrinsically disordered proteins. ``soursop.ped`` wraps PED's REST API so you can search it, look at entries, and load ensembles straight into SOURSOP as an :class:`~soursop.sstrajectory.SSTrajectory`, with one frame per conformer.

PED identifiers come in two forms. An **entry** (``PED00001``) groups one or more ensembles deposited together. An **ensemble** within that entry is written ``PED00001e001``. Every function that needs an ensemble accepts either the full form, or an entry plus ``ensemble_id`` (``"e001"`` or just ``1``). An entry on its own is accepted when it contains only one ensemble.

Quick start
----------------------------

.. code-block:: python

    from soursop import ped

    # search: returns entry identifiers, or ensemble identifiers
    ped.search("alpha-synuclein")
    ped.search(uniprot="P37840", ensembles=True)

    # an entry's metadata and the ensembles it contains
    entry = ped.get_entry("PED00001")
    print(entry)                 # PED00001: ... [3 ensemble(s): e001, e002, e003]
    print(entry.identifiers)     # ('PED00001e001', 'PED00001e002', 'PED00001e003')

    # load an ensemble and analyse it like any other trajectory
    traj = ped.load_ensemble("PED00001e001")
    protein = traj.proteinTrajectoryList[0]
    rg = protein.get_radius_of_gyration()

Searching
----------------------------

:func:`~soursop.ped.search` returns identifiers, and :func:`~soursop.ped.search_entries` returns full :class:`~soursop.ped.PEDEntry` records. Both take the same criteria, and every criterion you give must match:

* ``query`` - free text across the entry (e.g. ``"SAXS"``, ``"tau"``);
* ``protein_name``, ``uniprot`` - the protein;
* ``term`` - an IDPO term, i.e. a measurement or ensemble-generation method (e.g. ``"NMR"``);
* ``author``, ``publication`` (PubMed ID or DOI), ``publication_title``;
* ``cross_ref`` - a database or experimental cross-reference (e.g. a DisProt or BMRB identifier).

Results are paged through automatically. Use ``max_results`` to stop early. PED holds over 25,000 entries (mostly AlphaFlex and molecular dynamics ensembles), so a search with no criteria takes several minutes.

Two of PED's filters behave unexpectedly. ``author`` matches PED's *data owner* (the depositor) rather than the paper's authors, and ``cross_ref`` currently returns a server error from PED. A free-text search finds both, e.g. ``ped.search("Forman-Kay")`` or ``ped.search("DP00631")``.

Entries
----------------------------

:func:`~soursop.ped.get_entry` returns a :class:`~soursop.ped.PEDEntry` with the title, authors, UniProt accessions, protein names, methods and cross-references. It also holds one :class:`~soursop.ped.PEDEnsembleSummary` per ensemble, giving the number of conformers, the chain sequences, whether the ensemble is C-alpha only, and ``ped_statistics``. The last of these is PED's own analysis of the ensemble (mean Rg, end-to-end distance, asphericity, Flory exponent and so on), with distances in **nanometres**. SOURSOP's all-atom Rg reproduces PED's values exactly, which makes a handy check that an ensemble loaded correctly.

Loading ensembles
----------------------------

:func:`~soursop.ped.load_ensemble` downloads the ensemble's multi-model PDB file, reads it into an ``SSTrajectory`` and, by default, deletes the downloaded file. Pass ``cache_dir`` to keep the files, so each ensemble is only downloaded once. Any other keyword arguments go to ``SSTrajectory`` (e.g. ``protein_grouping``). To just save the PDB file, use :func:`~soursop.ped.download_ensemble`.

Some PED ensembles are **weighted**, meaning their conformers are not equally likely. :func:`~soursop.ped.get_ensemble_weights` returns the per-conformer weights, normalised to sum to 1, or ``None`` if the ensemble is unweighted. The result can be passed straight to SOURSOP's ``weights=`` arguments:

.. code-block:: python

    w = ped.get_ensemble_weights("PED00001e001")
    rg = protein.get_radius_of_gyration(weights=w)   # w=None -> unweighted

PED's own analyses
----------------------------

PED's API description lists per-ensemble analysis files (DSSP, radius of gyration and Dmax, Ramachandran and rotamer statistics, clash and geometry checks), and :func:`~soursop.ped.get_ensemble_asset` fetches any of them by name (see :data:`~soursop.ped.ENSEMBLE_ASSETS`). At the time of writing, however, PED no longer serves most of these as separate files: requests for the gyration and Dmax files return 404, and the DSSP file is empty. The summary statistics are available through ``ped_statistics`` instead. :func:`~soursop.ped.download_all_data` downloads everything PED holds for an ensemble (currently the ensemble and its weights) as a single ``.tar.gz``.

Errors and network access
----------------------------

Every request goes over HTTPS to ``https://deposition.proteinensemble.org/v1``. Failures raise :class:`~soursop.ped.PEDError`, a subclass of :class:`~soursop.ssexceptions.SSException`. Its ``status`` attribute holds the HTTP status code (e.g. ``404`` for an entry that does not exist), or ``None`` if PED could not be reached. Every function takes a ``timeout`` in seconds (default 120).

API reference
----------------------------

.. autofunction:: soursop.ped.search
.. autofunction:: soursop.ped.search_entries
.. autofunction:: soursop.ped.get_entry
.. autofunction:: soursop.ped.list_ensembles
.. autofunction:: soursop.ped.load_ensemble
.. autofunction:: soursop.ped.download_ensemble
.. autofunction:: soursop.ped.get_ensemble_weights
.. autofunction:: soursop.ped.get_ensemble_asset
.. autofunction:: soursop.ped.download_all_data
.. autofunction:: soursop.ped.parse_identifier
.. autofunction:: soursop.ped.full_identifier
.. autofunction:: soursop.ped.normalise_ensemble_id

.. autoclass:: soursop.ped.PEDEntry
   :members: from_json, ensemble_ids, identifiers, ensemble, load

.. autoclass:: soursop.ped.PEDEnsembleSummary
   :members: identifier, load

.. autoclass:: soursop.ped.PEDChainStats

.. autoclass:: soursop.ped.PEDError

.. autodata:: soursop.ped.ENSEMBLE_ASSETS
.. autodata:: soursop.ped.TABULAR_ASSETS
