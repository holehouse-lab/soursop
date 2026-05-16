"""Generate ground-truth reference dictionaries for the SOURSOP test suite.

For each entry in :data:`TRAJECTORIES`, this script loads the named trajectory,
computes every public observable on the ``SSProtein`` class that is meaningful
for the trajectory's resolution (AA = all-atom; CG = one-bead-per-residue),
and pickles the result to::

    soursop/data/test_data/<name>_<resolution>_reference.pkl

The pickles are meant to be consumed by a companion regression test that
recomputes the same observables on the same trajectory and compares each value
with ``numpy.testing.assert_allclose`` (or ``numpy.allclose``) to catch any
behavioural drift in SOURSOP.

Run this script ONLY when an observable's behaviour has changed deliberately
and the references should be refreshed. The output files should then be
committed together with the code change that motivated the refresh.

Usage
-----
    python soursop/tests/build_reference/build_references.py

Notes
-----
- Stochastic methods (``get_scaling_exponent`` and
  ``get_local_to_global_correlation``) are preceded by an explicit
  ``np.random.seed(SEED)`` so their output is reproducible. The full set of
  keys whose values were produced after a re-seed is recorded under
  ``ref['meta']['stochastic_keys']``.
- Per-frame observables are computed over the full trajectory but only
  ``SAVE_N_SNAPSHOTS`` (default: 10) evenly-spaced frames are saved. The save
  stride is therefore ``max(1, n_frames // SAVE_N_SNAPSHOTS)`` and is recorded
  in ``ref['meta']['save_stride']``. The future regression test must apply the
  same stride to its freshly computed outputs before comparison.
- Resolution-dependent skips: CG single-bead trajectories cannot produce
  meaningful values for dihedral angles, secondary structure (DSSP / BBSEG),
  sidechain alignment, SASA, native-contact Q, or any contact / inter-residue
  atomic-distance mode that depends on sidechain or heavy-atom selections. The
  corresponding keys are simply absent from CG reference dictionaries.
- Folded-protein caveat: scaling-exponent / polymer-scaled-distance-map /
  overlap-concentration assume polymer scaling. The numbers will be stable
  (good for regression testing) but are not physically meaningful for compact
  folded proteins.
"""

from __future__ import annotations

import math
import os
import pickle
import sys
from typing import Any, Callable

import numpy as np

import soursop
from soursop import __version__ as soursop_version
from soursop.sstrajectory import SSTrajectory
from soursop.ssexceptions import SSException
import mdtraj

# randseed for easy reproduciblity
SEED = 42

# Target number of per-frame snapshots to persist. The per-trajectory save
# stride is computed as ``max(1, n_frames // SAVE_N_SNAPSHOTS)``.
SAVE_N_SNAPSHOTS = 10

ANGLE_NAMES = ['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5']

# All contact-map modes that work for an all-atom trajectory. ``sidechain-heavy``
# is excluded because it raises explicitly when any GLY is present.
AA_CONTACT_MODES = ['closest-heavy', 'closest', 'ca', 'sidechain']

# For a single-bead-per-residue CG trajectory, ``ca`` is the only contact mode
# that is well-defined (every other mode depends on having more than one atom
# per residue or on a sidechain selection).
CG_CONTACT_MODES = ['ca']

POLYMER_SCALED_MAP_MODES = [
    'fractional-change',
    'signed-fractional-change',
    'signed-absolute-change',
    'scaled',
]

# Trajectories for which a reference pickle should be built. Each entry must
# specify the basename (which expands to ``<name>_<resolution>.{pdb,xtc}`` in
# ``soursop/data/test_data/``), the resolution tag (``AA`` or ``CG``), and the
# canonical R1/R2 residue indices used for inter-residue / regional analyses.
TRAJECTORIES = [
    {'name': 'ctl9', 'resolution': 'AA', 'R1': 5, 'R2': 45},
    {'name': 'sigA', 'resolution': 'CG', 'R1': 5, 'R2': 45},
    {'name': 'ntl9', 'resolution': 'AA', 'R1': 3, 'R2': 25},
    {'name': 'gs6', 'resolution': 'AA', 'R1': 2, 'R2': 4},
    {'name': 'synth_1', 'resolution': 'CG', 'R1': 2, 'R2': 40},
]


def _seed_then(func: Callable, *args, **kwargs):
    np.random.seed(SEED)
    return func(*args, **kwargs)


def _sub(arr, stride: int, axis: int = 0):
    """Subsample a per-frame array by ``stride`` along the given axis."""
    arr = np.asarray(arr)
    if arr.ndim == 0:
        return arr
    idx = np.arange(0, arr.shape[axis], stride)
    return np.take(arr, idx, axis=axis)


def _properties(p, resolution: str) -> dict:
    # CG trajectories without a periodic box have traj.unitcell_lengths == None;
    # the `unitcell` property raises TypeError in that case. Store a sentinel.
    try:
        unitcell = np.asarray(p.unitcell)
    except TypeError:
        unitcell = None
    return {
        'n_frames': int(p.n_frames),
        'n_residues': int(p.n_residues),
        'residue_index_list': list(p.residue_index_list),
        'resid_with_CA': list(p.resid_with_CA),
        'ncap': bool(p.ncap),
        'ccap': bool(p.ccap),
        'unitcell': unitcell,
        'sequence_three_letter': p.get_amino_acid_sequence(oneletter=False, numbered=False),
        'sequence_one_letter': p.get_amino_acid_sequence(oneletter=True, numbered=False),
        'resolution': resolution,
    }


def _global_structural(p, save_stride: int) -> dict:
    out: dict[str, Any] = {}
    out['get_radius_of_gyration'] = _sub(p.get_radius_of_gyration(), save_stride)
    out['get_center_of_mass'] = _sub(p.get_center_of_mass(), save_stride)
    out['get_gyration_tensor'] = _sub(p.get_gyration_tensor(verbose=False), save_stride)
    out['get_asphericity'] = _sub(p.get_asphericity(verbose=False), save_stride)
    out['get_hydrodynamic_radius__nygaard'] = _sub(
        p.get_hydrodynamic_radius(mode='nygaard'), save_stride,
    )
    out['get_hydrodynamic_radius__kr_CA'] = _sub(
        p.get_hydrodynamic_radius(mode='kr', distance_mode='CA'), save_stride,
    )
    out['get_hydrodynamic_radius__kr_COM'] = _sub(
        p.get_hydrodynamic_radius(mode='kr', distance_mode='COM'), save_stride,
    )
    out['get_t'] = _sub(p.get_t(), save_stride)
    out['get_molecular_volume'] = _sub(p.get_molecular_volume(), save_stride)
    out['get_end_to_end_distance__CA'] = _sub(p.get_end_to_end_distance(mode='CA'), save_stride)
    out['get_end_to_end_distance__COM'] = _sub(p.get_end_to_end_distance(mode='COM'), save_stride)
    out['get_overlap_concentration'] = float(p.get_overlap_concentration())
    out['get_end_to_end_vs_rg_correlation__CA'] = float(
        p.get_end_to_end_vs_rg_correlation(mode='CA')
    )
    out['get_end_to_end_vs_rg_correlation__COM'] = float(
        p.get_end_to_end_vs_rg_correlation(mode='COM')
    )
    return out


def _distance(p, R1: int, R2: int, save_stride: int, resolution: str) -> dict:
    out: dict[str, Any] = {}
    out['calculate_all_CA_distances__CA'] = _sub(
        p.calculate_all_CA_distances(R1, mode='CA'), save_stride,
    )
    out['calculate_all_CA_distances__COM'] = _sub(
        p.calculate_all_CA_distances(R1, mode='COM'), save_stride,
    )

    # distance maps are aggregated over all frames (mean/std) — keep full
    mean_CA, std_CA = p.get_distance_map(mode='CA', verbose=False)
    out['get_distance_map__CA'] = {'mean': mean_CA, 'std': std_CA}

    mean_COM, std_COM = p.get_distance_map(mode='COM', verbose=False)
    out['get_distance_map__COM'] = {'mean': mean_COM, 'std': std_COM}

    mean_CA_RMS, std_CA_RMS = p.get_distance_map(mode='CA', RMS=True, verbose=False)
    out['get_distance_map__CA_RMS'] = {'mean': mean_CA_RMS, 'std': std_CA_RMS}

    out['get_inter_residue_COM_distance'] = _sub(
        p.get_inter_residue_COM_distance(R1, R2), save_stride,
    )
    out['get_inter_residue_COM_vector'] = _sub(
        p.get_inter_residue_COM_vector(R1, R2), save_stride,
    )
    out['get_inter_residue_atomic_distance__atom_CA_CA'] = _sub(
        p.get_inter_residue_atomic_distance(R1, R2, A1='CA', A2='CA', mode='atom'),
        save_stride,
    )
    out['get_inter_residue_atomic_distance__ca'] = _sub(
        p.get_inter_residue_atomic_distance(R1, R2, mode='ca'), save_stride,
    )

    if resolution == 'AA':
        # closest / closest-heavy / sidechain modes all require multiple atoms
        # per residue, which CG single-bead topologies do not provide.
        out['get_inter_residue_atomic_distance__closest'] = _sub(
            p.get_inter_residue_atomic_distance(R1, R2, mode='closest'), save_stride,
        )
        out['get_inter_residue_atomic_distance__closest_heavy'] = _sub(
            p.get_inter_residue_atomic_distance(R1, R2, mode='closest-heavy'), save_stride,
        )
        out['get_inter_residue_atomic_distance__sidechain'] = _sub(
            p.get_inter_residue_atomic_distance(R1, R2, mode='sidechain'), save_stride,
        )
    return out


def _angles_dihedrals(p, save_stride: int) -> dict:
    """AA-only: dihedral angles require backbone N/C/CA and sidechain atoms."""
    out: dict[str, Any] = {}

    for name in ANGLE_NAMES:
        atom_names, angles = p.get_angles(name)
        # atom_names are mdtraj Atom objects; stringify so the pickle is portable
        stringified = [[str(a) for a in residue_atoms] for residue_atoms in atom_names]
        # angles shape is (nres, n_frames) — subsample along the frame axis
        out[f'get_angles__{name}'] = {
            'atom_names': stringified,
            'angles': _sub(np.asarray(angles), save_stride, axis=1),
        }

    out['get_angle_decay__C_N'] = p.get_angle_decay(atom1='C', atom2='N')

    out['get_dihedral_mutual_information__psi__norm0'] = p.get_dihedral_mutual_information(
        angle_name='psi', normalize=False,
    )
    out['get_dihedral_mutual_information__psi__norm1'] = p.get_dihedral_mutual_information(
        angle_name='psi', normalize=True,
    )
    return out


def _polymer(p, stochastic_keys: list) -> dict:
    out: dict[str, Any] = {}

    # mean_vals=True keeps the per-separation mean only; the raw per-frame
    # distances add no extra regression coverage beyond what
    # get_internal_scaling_RMS already provides.
    seq_CA, mean_CA = p.get_internal_scaling(mode='CA', mean_vals=True, verbose=False)
    out['get_internal_scaling__CA'] = {
        'seq_sep': np.asarray(seq_CA),
        'mean': np.asarray(mean_CA),
    }

    seq_COM, mean_COM = p.get_internal_scaling(mode='COM', mean_vals=True, verbose=False)
    out['get_internal_scaling__COM'] = {
        'seq_sep': np.asarray(seq_COM),
        'mean': np.asarray(mean_COM),
    }

    seq_CA_rms, mean_CA_rms = p.get_internal_scaling_RMS(mode='CA', verbose=False)
    out['get_internal_scaling_RMS__CA'] = {
        'seq_sep': np.asarray(seq_CA_rms),
        'mean': np.asarray(mean_CA_rms),
    }

    seq_COM_rms, mean_COM_rms = p.get_internal_scaling_RMS(mode='COM', verbose=False)
    out['get_internal_scaling_RMS__COM'] = {
        'seq_sep': np.asarray(seq_COM_rms),
        'mean': np.asarray(mean_COM_rms),
    }

    se_keys = ['best_nu', 'best_A0', 'min_nu', 'max_nu', 'min_A0', 'max_A0',
               'redchi_fit_region', 'redchi_all_points',
               'fit_region_data', 'all_points_data']

    # get_scaling_exponent and get_polymer_scaled_distance_map fit a homopolymer
    # model and need enough residues to do so (at least ~25). For very short
    # chains (e.g., gs6 with 6 residues) the fit isn't well-defined and the
    # function raises. Skip these observables gracefully so the rest of the
    # reference still builds; the absent keys are then skipped by smoke_check.
    try:
        se_CA = _seed_then(p.get_scaling_exponent, mode='CA', verbose=False)
        out['get_scaling_exponent__CA'] = dict(zip(se_keys, se_CA))
        stochastic_keys.append('polymer.get_scaling_exponent__CA')

        se_COM = _seed_then(p.get_scaling_exponent, mode='COM', verbose=False)
        out['get_scaling_exponent__COM'] = dict(zip(se_keys, se_COM))
        stochastic_keys.append('polymer.get_scaling_exponent__COM')

        for mode in POLYMER_SCALED_MAP_MODES:
            psd_map, nu, A0, redchi = _seed_then(
                p.get_polymer_scaled_distance_map, mode=mode, verbose=False,
            )
            out[f'get_polymer_scaled_distance_map__{mode}'] = {
                'map': psd_map,
                'nu': float(nu),
                'A0': float(A0),
                'redchi': float(redchi),
            }
            stochastic_keys.append(f'polymer.get_polymer_scaled_distance_map__{mode}')
    except SSException as e:
        print(
            f"  [WARNING] skipping polymer-scaling fits (chain too short for "
            f"homopolymer model): {e}"
        )

    return out


def _secondary_structure(p, save_stride: int) -> dict:
    """AA-only: DSSP needs full backbone, BBSEG needs phi/psi."""
    out: dict[str, Any] = {}

    dssp_sum = p.get_secondary_structure_DSSP(return_per_frame=False)
    out['DSSP_summary'] = {
        'resid_list': np.asarray(dssp_sum[0]),
        'helix': dssp_sum[1],
        'extended': dssp_sum[2],
        'coil': dssp_sum[3],
    }

    dssp_pf = p.get_secondary_structure_DSSP(return_per_frame=True)
    out['DSSP_per_frame'] = {
        'resid_list': np.asarray(dssp_pf[0]),
        'helix': _sub(dssp_pf[1], save_stride),
        'extended': _sub(dssp_pf[2], save_stride),
        'coil': _sub(dssp_pf[3], save_stride),
    }

    bbseg_sum = p.get_secondary_structure_BBSEG(return_per_frame=False)
    out['BBSEG_summary'] = {
        'resid_list': np.asarray(bbseg_sum[0]),
        'per_class': {int(k): np.asarray(v) for k, v in bbseg_sum[1].items()},
    }

    bbseg_pf = p.get_secondary_structure_BBSEG(return_per_frame=True)
    out['BBSEG_per_frame'] = {
        'resid_list': np.asarray(bbseg_pf[0]),
        'per_class': {int(k): _sub(np.asarray(v), save_stride) for k, v in bbseg_pf[1].items()},
    }
    return out


def _rmsd(p, save_stride: int) -> dict:
    return {
        'get_RMSD__frame0_to_lastframe': p.get_RMSD(frame1=0, frame2=p.n_frames - 1),
        'get_RMSD__frame0_to_all_frames': _sub(p.get_RMSD(frame1=0, frame2=-1), save_stride),
    }


def _contacts(p, save_stride: int, resolution: str) -> dict:
    out: dict[str, Any] = {}

    modes = AA_CONTACT_MODES if resolution == 'AA' else CG_CONTACT_MODES
    for mode in modes:
        cmap, corder = p.get_contact_map(mode=mode)
        key = mode.replace('-', '_')
        out[f'get_contact_map__{key}'] = {
            'contact_map': cmap,
            'contact_order': corder,
        }

    if resolution != 'AA':
        # get_Q relies on heavy-atom contacts which are undefined for CG.
        return out

    # get_Q requires at least one native contact (heavy-atom pair within
    # threshold and separated by >3 residues in sequence). For very short
    # peptides this can be zero, in which case get_Q returns NaN. Skip
    # gracefully in that case.
    q_avg = p.get_Q(protein_average=True)
    if not np.all(np.isfinite(np.asarray(q_avg))):
        print(
            "  [WARNING] skipping get_Q: native-contact count is zero "
            "(chain too short or too extended for any native contacts)"
        )
        return out

    out['get_Q__protein_average_True'] = _sub(q_avg, save_stride)

    # get_Q(protein_average=False) returns a 5-tuple:
    #   [0] per-native-contact fractional native (vector)
    #   [1] native contact atom-pair definitions (n x 2 int array)
    #   [2] residue-keyed dict {resname-resid: list[float]}
    #   [3] ordered keys from [2]
    #   [4] nres x nres inter-residue Q matrix
    q_per_contact, native_pairs, res2contacts, ordered_keys, res_res_q = p.get_Q(
        protein_average=False,
    )
    out['get_Q__protein_average_False'] = {
        'per_contact_fraction': np.asarray(q_per_contact),
        'native_contact_pairs': np.asarray(native_pairs),
        'per_residue_contacts': {k: np.asarray(v) for k, v in res2contacts.items()},
        'ordered_residue_keys': list(ordered_keys),
        'res_res_q_matrix': np.asarray(res_res_q),
    }
    return out


def _sasa(p, R1: int, R2: int, save_stride: int) -> dict:
    """AA-only: a single CA atom per residue is not a meaningful SASA."""
    out: dict[str, Any] = {}
    out['get_all_SASA__residue'] = p.get_all_SASA(mode='residue', stride=save_stride)
    out['get_all_SASA__atom'] = p.get_all_SASA(mode='atom', stride=save_stride)
    out['get_regional_SASA'] = p.get_regional_SASA(R1, R2, stride=save_stride)

    site = p.get_site_accessibility(
        ['ALA', 'LEU', 'VAL'], mode='residue_type', stride=save_stride,
    )
    out['get_site_accessibility__ALA_LEU_VAL'] = {
        str(k): np.asarray(v) for k, v in site.items()
    }
    return out


def _alignment_heterogeneity(
    p, R1: int, R2: int, save_stride: int, resolution: str, stochastic_keys: list,
) -> dict:
    out: dict[str, Any] = {}

    if resolution == 'AA':
        out['get_sidechain_alignment_angle'] = _sub(
            p.get_sidechain_alignment_angle(R1, R2), save_stride,
        )

    out['get_D_vector'] = p.get_D_vector(stride=save_stride, verbose=False)

    # get_local_heterogeneity / get_local_collapse use a sliding-window over
    # the chain and require n_residues >= fragment_size/window_size. Skip
    # gracefully if the chain is shorter than the default 10 residues.
    try:
        lh_mean, lh_std, lh_histo, lh_bins = p.get_local_heterogeneity(
            fragment_size=10, stride=save_stride, verbose=False,
        )
        out['get_local_heterogeneity'] = {
            'mean': np.asarray(lh_mean),
            'std': np.asarray(lh_std),
            'histo': np.asarray(lh_histo),
            'bins': np.asarray(lh_bins),
        }
    except SSException as e:
        print(f"  [WARNING] skipping get_local_heterogeneity: {e}")

    try:
        lc_mean, lc_std, lc_histo, lc_bins = p.get_local_collapse(
            window_size=10, verbose=False,
        )
        out['get_local_collapse'] = {
            'mean': np.asarray(lc_mean),
            'std': np.asarray(lc_std),
            'histo': np.asarray(lc_histo),
            'bins': np.asarray(lc_bins),
        }
    except SSException as e:
        print(f"  [WARNING] skipping get_local_collapse: {e}")

    # get_local_to_global_correlation needs enough residues to vary sequence
    # separation. Skip gracefully (along with its stochastic_keys entry) if
    # it raises.
    try:
        raw, n_pairs, mean_corr, std_corr = _seed_then(
            p.get_local_to_global_correlation,
            mode='COM', n_cycles=100, max_num_pairs=10, stride=save_stride, verbose=False,
        )
        out['get_local_to_global_correlation'] = {
            'raw': raw,
            'n_pairs': n_pairs,
            'mean': mean_corr,
            'std': std_corr,
        }
        stochastic_keys.append('alignment_heterogeneity.get_local_to_global_correlation')
    except SSException as e:
        print(f"  [WARNING] skipping get_local_to_global_correlation: {e}")

    return out


def _validate_canonical_residues(p, R1: int, R2: int, resolution: str) -> None:
    seq = p.get_amino_acid_sequence(oneletter=False, numbered=False)
    # Caps must never be selected. GLY must be avoided for AA because
    # get_sidechain_alignment_angle raises on glycine; for CG that method is
    # skipped, so glycines are tolerated.
    forbidden = {'ACE', 'NME'} | ({'GLY'} if resolution == 'AA' else set())
    for label, resid in (('R1', R1), ('R2', R2)):
        if resid >= len(seq):
            raise RuntimeError(
                f"{label}={resid} is out of bounds for sequence length {len(seq)}"
            )
        resname = seq[resid].split('-')[0]
        if resname in forbidden:
            raise RuntimeError(
                f"{label}={resid} is {resname}; not allowed for resolution={resolution}. "
                "Update R1/R2 in the TRAJECTORIES entry."
            )


def build_reference(protein, name: str, resolution: str, R1: int, R2: int) -> dict:
    """Compute every observable and assemble the reference dictionary."""
    _validate_canonical_residues(protein, R1, R2, resolution)
    save_stride = max(1, int(protein.n_frames) // SAVE_N_SNAPSHOTS)

    stochastic_keys: list[str] = []

    properties = _properties(protein, resolution)
    global_structural = _global_structural(protein, save_stride)
    distance = _distance(protein, R1, R2, save_stride, resolution)
    polymer = _polymer(protein, stochastic_keys)
    rmsd = _rmsd(protein, save_stride)
    contacts = _contacts(protein, save_stride, resolution)
    alignment = _alignment_heterogeneity(
        protein, R1, R2, save_stride, resolution, stochastic_keys,
    )

    meta = {
        'name': name,
        'resolution': resolution,
        'soursop_version': soursop_version,
        'mdtraj_version': mdtraj.__version__,
        'numpy_version': np.__version__,
        'python_version': '.'.join(str(v) for v in sys.version_info[:3]),
        'n_frames': properties['n_frames'],
        'n_residues': properties['n_residues'],
        'seed': SEED,
        'R1': R1,
        'R2': R2,
        'save_stride': save_stride,
        'save_n_snapshots_target': SAVE_N_SNAPSHOTS,
        'stochastic_keys': sorted(stochastic_keys),
    }

    ref: dict[str, Any] = {
        'meta': meta,
        'properties': properties,
        'global': global_structural,
        'distance': distance,
        'polymer': polymer,
        'rmsd': rmsd,
        'contacts': contacts,
        'alignment_heterogeneity': alignment,
    }

    if resolution == 'AA':
        ref['angles'] = _angles_dihedrals(protein, save_stride)
        ref['secondary_structure'] = _secondary_structure(protein, save_stride)
        ref['sasa'] = _sasa(protein, R1, R2, save_stride)

    return ref


def _walk_leaves(d: dict, prefix: str = ''):
    for key, value in d.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            yield from _walk_leaves(value, path)
        else:
            yield path, value


NAN_INF_WHITELIST: set[str] = set()


def _array_summary(value) -> str:
    if isinstance(value, np.ndarray):
        return f"shape={value.shape} dtype={value.dtype}"
    if isinstance(value, (list, tuple)):
        return f"len={len(value)} type={type(value).__name__}"
    if isinstance(value, (int, float, np.floating, np.integer)):
        return f"scalar={value!r}"
    if isinstance(value, bool):
        return f"bool={value}"
    if isinstance(value, str):
        return f"str(len={len(value)})"
    return f"type={type(value).__name__}"


def smoke_check(ref: dict, protein) -> None:
    print('--- smoke_check: walking reference dict ---')
    leaves = list(_walk_leaves(ref))
    for path, value in leaves:
        print(f"  {path:65s} {_array_summary(value)}")

    print(f'--- smoke_check: {len(leaves)} leaf values ---')

    # 'properties.unitcell' is legitimately None for CG trajectories without
    # a periodic box; any other None leaf is a bug.
    for path, value in leaves:
        if value is None and path != 'properties.unitcell':
            raise AssertionError(f"leaf '{path}' is None")

    for path, value in leaves:
        if path in NAN_INF_WHITELIST:
            continue
        if isinstance(value, np.ndarray) and value.dtype.kind in ('f', 'c'):
            if not np.all(np.isfinite(value)):
                bad = np.sum(~np.isfinite(value))
                raise AssertionError(
                    f"leaf '{path}' contains {bad} non-finite (NaN/Inf) entries"
                )
        elif isinstance(value, float) and not math.isfinite(value):
            raise AssertionError(f"leaf '{path}' is non-finite scalar: {value}")

    # Stochastic keys vary by chain length: very short chains (e.g., gs6 with
    # 6 CA-bearing residues) skip both polymer-scaling fits AND
    # get_local_to_global_correlation, ending up with zero stochastic keys.
    # That's a legitimate outcome, not a bug, so we no longer require >= 1.

    print('--- smoke_check: re-running stochastic methods to verify seeding ---')

    se_keys = ['best_nu', 'best_A0', 'min_nu', 'max_nu', 'min_A0', 'max_A0',
               'redchi_fit_region', 'redchi_all_points',
               'fit_region_data', 'all_points_data']

    if 'get_scaling_exponent__CA' in ref['polymer']:
        se_CA_2 = _seed_then(protein.get_scaling_exponent, mode='CA', verbose=False)
        se_CA_2_dict = dict(zip(se_keys, se_CA_2))
        for k, v in ref['polymer']['get_scaling_exponent__CA'].items():
            np.testing.assert_allclose(
                np.asarray(se_CA_2_dict[k]), np.asarray(v),
                err_msg=f"stochastic re-run mismatch: get_scaling_exponent__CA.{k}",
            )

    if 'get_scaling_exponent__COM' in ref['polymer']:
        se_COM_2 = _seed_then(protein.get_scaling_exponent, mode='COM', verbose=False)
        se_COM_2_dict = dict(zip(se_keys, se_COM_2))
        for k, v in ref['polymer']['get_scaling_exponent__COM'].items():
            np.testing.assert_allclose(
                np.asarray(se_COM_2_dict[k]), np.asarray(v),
                err_msg=f"stochastic re-run mismatch: get_scaling_exponent__COM.{k}",
            )

    for mode in POLYMER_SCALED_MAP_MODES:
        key = f'get_polymer_scaled_distance_map__{mode}'
        if key not in ref['polymer']:
            continue
        psd_2 = _seed_then(
            protein.get_polymer_scaled_distance_map, mode=mode, verbose=False,
        )
        ref_entry = ref['polymer'][key]
        np.testing.assert_allclose(psd_2[0], ref_entry['map'])
        assert float(psd_2[1]) == ref_entry['nu']
        assert float(psd_2[2]) == ref_entry['A0']
        assert float(psd_2[3]) == ref_entry['redchi']

    if 'get_local_to_global_correlation' in ref['alignment_heterogeneity']:
        # Must use the same stride the build call used (stored in meta);
        # otherwise the random-number consumption pattern diverges and the
        # re-run won't match.
        save_stride = ref['meta']['save_stride']
        raw2, n_pairs2, mean2, std2 = _seed_then(
            protein.get_local_to_global_correlation,
            mode='COM', n_cycles=100, max_num_pairs=10, stride=save_stride, verbose=False,
        )
        ref_l2g = ref['alignment_heterogeneity']['get_local_to_global_correlation']
        np.testing.assert_allclose(raw2, ref_l2g['raw'])
        np.testing.assert_allclose(n_pairs2, ref_l2g['n_pairs'])
        np.testing.assert_allclose(mean2, ref_l2g['mean'])
        np.testing.assert_allclose(std2, ref_l2g['std'])

    print('--- smoke_check: all stochastic re-runs match ---')


def _build_one(name: str, resolution: str, R1: int, R2: int) -> None:
    pdb = soursop.get_data(f'test_data/{name}_{resolution}.pdb')
    xtc = soursop.get_data(f'test_data/{name}_{resolution}.xtc')
    output_path = soursop.get_data(f'test_data/{name}_{resolution}_reference.pkl')

    print(f"\n===== {name}_{resolution} =====")
    print(f"Loading trajectory: pdb={pdb}\n                    xtc={xtc}")
    traj = SSTrajectory(xtc, pdb)
    protein = traj.proteinTrajectoryList[0]
    print(f"Loaded SSProtein with {protein.n_residues} residues, {protein.n_frames} frames")

    ref = build_reference(protein, name, resolution, R1, R2)
    smoke_check(ref, protein)

    print(f"Writing reference pickle to: {output_path}")
    with open(output_path, 'wb') as f:
        pickle.dump(ref, f, protocol=pickle.HIGHEST_PROTOCOL)

    size_bytes = os.path.getsize(output_path)
    print(f"Wrote {size_bytes:,} bytes ({size_bytes / (1024 * 1024):.2f} MiB)")


def main() -> None:
    for entry in TRAJECTORIES:
        _build_one(
            name=entry['name'],
            resolution=entry['resolution'],
            R1=entry['R1'],
            R2=entry['R2'],
        )
    print("\nAll references built.")


if __name__ == '__main__':
    main()
