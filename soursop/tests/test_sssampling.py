"""
Unit and regression test for the soursop package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
from soursop.sssampling import hellinger_distance, SamplingQuality
from soursop.sstools import find_trajectory_files
import soursop
import os


def test_hellingers_distance():
    p = np.array([1,0])
    q = np.array([0,1])
    
    assert hellinger_distance(p,q) == 1.0 

    p = np.array([0,0])
    q = np.array([0,0])
    assert hellinger_distance(p,q) == 0.0 

    p = np.array([[1,0],[1,0]])
    q = np.array([[0,1],[0,1]])
    assert np.all(hellinger_distance(p,q) == np.array([1.,1.]))

    p = np.array([[0,0],[1,0]])
    q = np.array([[0,0],[0,1]])
    assert np.all(hellinger_distance(p,q) == np.array([0.,1.]))

def test_sampling_quality():
    data_dir = soursop.get_data('test_data')
    wt_data = os.path.join(data_dir,"sampling_quality/WT")
    ev_data = os.path.join(data_dir,"sampling_quality/EV")
    wt_traj_paths, wt_top_paths  = find_trajectory_files(wt_data,3)
    ev_traj_paths, ev_top_paths = find_trajectory_files(ev_data,3)
    
    alanine_qual = SamplingQuality(wt_traj_paths,
                    ev_traj_paths, 
                    top_file=wt_top_paths,
                    ref_top=ev_top_paths,
                    method="1D angle distributions")
    
    hellingers = np.array([[[0.24167835, 0.19611402],
                            [0.24688875, 0.24731097],
                            [0.24903448, 0.20723213]],
                            [[0.13862987, 0.1848602 ],
                             [0.159583  , 0.18751601],
                             [0.1429838 , 0.13090839]]])
    
    assert np.all(alanine_qual.psi_angles.shape == (3, 2, 500))
    assert np.allclose(alanine_qual.psi_angles[0][0][0], -39.17526)
    assert np.allclose(alanine_qual.phi_angles[0][0][0], -144.49619)
    assert np.allclose(alanine_qual.compute_dihedral_hellingers(), hellingers)

    bins = alanine_qual.get_degree_bins()
    test_traj = alanine_qual.trajs[0]
    test_ref_traj = alanine_qual.ref_trajs[0]

    test_psi = test_traj.proteinTrajectoryList[0].get_angles("psi")[1]
    test_ref_psi = test_ref_traj.proteinTrajectoryList[0].get_angles("psi")[1]

    test_pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0],
                                     axis=1, 
                                     arr=test_psi)*np.round(np.rad2deg(alanine_qual.bwidth))

    test_ref_pdf = np.apply_along_axis(lambda col: np.histogram(col, bins=bins, density=True)[0],
                                         axis=1,
                                         arr=test_ref_psi)*np.round(np.rad2deg(alanine_qual.bwidth))

    assert np.allclose(
                        test_pdf,
                        alanine_qual.compute_pdf(alanine_qual.psi_angles[0],
                        bins=bins)
                    )

    assert np.allclose(
                        test_ref_pdf, 
                        alanine_qual.compute_pdf(alanine_qual.ref_psi_angles[0], 
                        bins=bins)
                    )



def test_hellingers_2d():
    data_dir = soursop.get_data('test_data')
    wt_data = os.path.join(data_dir,"sampling_quality/WT")
    ev_data = os.path.join(data_dir,"sampling_quality/EV")
    wt_traj_paths, wt_top_paths  = find_trajectory_files(wt_data,3)
    ev_traj_paths, ev_top_paths = find_trajectory_files(ev_data,3)
    
    hellingers_2d = np.array([[0.44302328, 0.42337237],
                            [0.46023975, 0.46853325],
                            [0.42439771, 0.42886934]])
    
    alanine_qual = SamplingQuality(wt_traj_paths,
                    ev_traj_paths, 
                    top_file=wt_top_paths,
                    ref_top=ev_top_paths,
                    method="2D angle distributions")
    
    assert np.allclose(alanine_qual.compute_dihedral_hellingers(), hellingers_2d)