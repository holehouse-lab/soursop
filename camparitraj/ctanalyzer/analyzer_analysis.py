##
##                                       _ _              _
##   ___ __ _ _ __ ___  _ __   __ _ _ __(_) |_ _ __ __ _ (_)
##  / __/ _` | '_ ` _ \| '_ \ / _` | '__| | __| '__/ _` || |
## | (_| (_| | | | | | | |_) | (_| | |  | | |_| | | (_| || |
##  \___\__,_|_| |_| |_| .__/ \__,_|_|  |_|\__|_|  \__,_|/ |
##                     |_|                             |__/
##
## Alex Holehouse (Pappu Lab and Holehouse Lab)
## Simulation analysis package
## Copyright 2014 - 2020
##

import numpy as np
import mdtraj as md
from . import configs
from camparitraj._internal_data import MAX_SASA_DATA
from .analyzer_exception import AnalyzerException


AALIST = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


def status_message(ana, name):
    print("Running %s on %s" % (ana, name))


def arrayfy(val):
    return np.array([val])


def run_RG(CP, outdir):
    status_message('Radius of gyration', outdir)
    RG = CP.get_radius_of_gyration()
    MEAN_RG = arrayfy(np.mean(RG))
    STD_RG  = arrayfy(np.std(RG))
    np.savetxt('%s/RG.csv' % (outdir), RG, delimiter=',')
    np.savetxt('%s/RG_mean.csv' % (outdir), MEAN_RG, delimiter=', ')
    np.savetxt('%s/RG_std.csv' % (outdir), STD_RG, delimiter=', ')


def run_t_inst(CP, outdir):
    status_message('Instantaneous t', outdir)
    T = CP.get_t()
    MEAN_T = arrayfy(np.mean(T))
    STD_T  = arrayfy(np.std(T))
    np.savetxt('%s/t_inst.csv' % (outdir), T, delimiter=',')
    np.savetxt('%s/t_mean.csv' % (outdir), MEAN_T, delimiter=', ')
    np.savetxt('%s/t_std.csv' % (outdir), STD_T, delimiter=', ')


def run_RH(CP, outdir):
    status_message('Hydrodynamic radius', outdir)
    RH = CP.get_hydrodynamic_radius()
    MEAN_RH = arrayfy(np.mean(RH))
    STD_RH  = arrayfy(np.std(RH))
    np.savetxt('%s/RH.csv' % (outdir), RH, delimiter=',')
    np.savetxt('%s/RH_mean.csv' % (outdir), MEAN_RH, delimiter=', ')
    np.savetxt('%s/RH_std.csv' % (outdir), STD_RH, delimiter=', ')


def run_end_to_end(CP, outdir):
    status_message('End to end distance', outdir)
    E2E = CP.get_end_to_end_distance()
    MEAN_E2E = arrayfy(np.mean(E2E))
    STD_E2E  = arrayfy(np.std(E2E))
    np.savetxt('%s/end_to_end.csv' % (outdir), E2E, delimiter=',')
    np.savetxt('%s/end_to_end_mean.csv' % (outdir), MEAN_E2E, delimiter=', ')
    np.savetxt('%s/end_to_end_std.csv' % (outdir), STD_E2E, delimiter=', ')


def run_asphericity(CP, outdir):
    status_message('Asphericity', outdir)
    asph = CP.get_asphericity()
    MEAN_asph = arrayfy(np.mean(asph))
    STD_asph  = arrayfy(np.std(asph))
    np.savetxt('%s/ASPH.csv' % (outdir), asph, delimiter=',')
    np.savetxt('%s/ASPH_mean.csv' % (outdir), MEAN_asph, delimiter=', ')
    np.savetxt('%s/ASPH_std.csv' % (outdir), STD_asph, delimiter=', ')


def run_distanceMap(CP, outdir):
    status_message('Distance map', outdir)
    [a,b] = CP.get_distance_map()
    np.savetxt('%s/distance_map.csv' % (outdir), a, delimiter=',')
    np.savetxt('%s/distance_map_std.csv' % (outdir), b, delimiter=', ')


def run_polymer_scaling_map(CP, outdir):
    status_message('Polymer scaling map', outdir)
    print("... absolute change:")
    [PSM, nu, A0, redchi] = CP.get_polymer_scaled_distance_map(mode='signed-absolute-change')
    np.savetxt('%s/polymer_deviation_map_absolute.csv' % (outdir), PSM, delimiter=',')

    print("...fractional change:")
    [PSM, nu, A0, redchi] = CP.get_polymer_scaled_distance_map(mode='signed-fractional-change')
    np.savetxt('%s/polymer_deviation_map_fractional.csv' % (outdir), PSM, delimiter=',')

    np.savetxt('%s/polymer_deviation_map_params.csv' % (outdir), np.transpose([nu,A0,redchi]), delimiter=', ')


def run_analytical_frc(CP, outdir,count=False):
    # This seems to be a missing module / package? It seems to be user-derived
    # as there is no record of it on the visible web.
    import afrc
    AAS = CP.get_amino_acid_sequence(oneletter=True)
    AAS_final = AAS.translate(str.maketrans('','','<>'))
    AFRC = afrc.AnalyticalFRC(AAS_final)

    if count == False:
        count = 50000

    # distance map


    # end_to_end
    [a,b] = AFRC.get_re_distribution()
    np.savetxt('%s/AFRC_end_to_end_distribution.csv' % (outdir),np.array((a,b)).transpose())

    re = AFRC.sample_re_distribution(n=count)
    np.savetxt('%s/AFRC_end_to_end.csv' % (outdir),re)

    [a,b] = AFRC.get_rg_distribution()
    np.savetxt('%s/AFRC_rg_distribution.csv' % (outdir),np.array((a,b)).transpose())

    rg = AFRC.sample_rg_distribution(n=count)
    np.savetxt('%s/AFRC_RG.csv' % (outdir),rg)

    dm = AFRC.get_distance_map()
    np.savetxt('%s/AFRC_distance_map.csv' % (outdir),dm)


def run_internal_scaling(CP, outdir):
    IS = CP.get_internal_scaling()
    mean_is = [np.mean(i) for i in IS[1]]
    np.savetxt('%s/INTSCAL.csv' % (outdir), mean_is, delimiter=', ')


def run_contact_map(CP, outdir, d_thresh):
    cmap_full = CP.get_contact_map(distance_thresh=d_thresh)
    np.savetxt('%s/contact_map_%3.3f.csv' % (outdir, d_thresh),cmap_full[0])
    np.savetxt('%s/contact_order_%3.3f.csv' % (outdir, d_thresh),cmap_full[1])


def run_RMS_internal_scaling(CP, outdir):
    IS = CP.get_internal_scaling_RMS()
    np.savetxt('%s/RMS_INTSCAL.csv' % (outdir), IS[1], delimiter=', ')


def run_fractal_deviation(CP, outdir, stride):
    (_, n_pairs, mean_cor, std_cor) = CP.get_local_to_global_correlation(stride=stride, n_cycles=1000, max_num_pairs=20)
    np.savetxt('%s/fractional_deviation.csv' % (outdir), np.transpose([n_pairs, mean_cor, std_cor]), delimiter=', ')


def run_Q_analysis(CP, outdir):
    Q_TUPLE =  CP.get_Q(stride=1, protein_average=False)
    with open('%s/Q_res_by_res.txt' % (outdir), 'w') as fh:
        for i in Q_TUPLE[3]:
            fh.write("%s, %3.3f\n" % (i[3:], np.mean(Q_TUPLE[2][i])))


def run_rij_analysis(CP, outdir, ri, rj):
    rij = CP.get_inter_residue_COM_distance(ri, rj)

    MEAN_rij = arrayfy(np.mean(rij))
    STD_rij  = arrayfy(np.std(rij))
    np.savetxt('%s/r_%i_%i.csv' % (outdir, ri, rj), rij, delimiter=',')
    np.savetxt('%s/r_%i_%i_mean.csv' % (outdir, ri, rj), MEAN_rij, delimiter=', ')
    np.savetxt('%s/r_%i_%i_std.csv' % (outdir, ri, rj), STD_rij, delimiter=', ')


def run_rg_re_correlation(CP, outdir):
    c = CP.get_end_to_end_vs_rg_correlation()
    np.savetxt('%s/rg_re_corr.csv' % (outdir), arrayfy(c), delimiter=', ')


def run_scaling_exponent_power(CP, outdir, end_effect=5):
    c = CP.get_scaling_exponent(end_effect=end_effect, mode='COM')

    if end_effect == configs.DEFAULT_END_EFFECT:
        outname_1='scaling_exp_analysis_power.csv'
        outname_2='scaling_exp_idx_used_power.csv'
        outname_3='scaling_exp_fit_power.csv'
    else:
        outname_1='scaling_exp_analysis_ee%i_power.csv' %(end_effect)
        outname_2='scaling_exp_idx_used_ee%i_power.csv' % (end_effect)
        outname_3='scaling_exp_fit_ee%i_power.csv' % (end_effect)

    np.savetxt('%s/%s' % (outdir, outname_1), c[0:8], delimiter=', ')
    np.savetxt('%s/%s' % (outdir, outname_2), c[8].transpose(), delimiter=', ')
    np.savetxt('%s/%s' % (outdir, outname_3), c[9].transpose(), delimiter=', ')


def run_scaling_exponent_power_CA(CP, outdir, end_effect=5):
    c = CP.get_scaling_exponent(end_effect=end_effect, mode='CA')

    if end_effect == configs.DEFAULT_END_EFFECT:
        outname_1='scaling_exp_analysis_power_CA.csv'
        outname_2='scaling_exp_idx_used_power_CA.csv'
        outname_3='scaling_exp_fit_power_CA.csv'
    else:
        outname_1='scaling_exp_analysis_ee%i_power_CA.csv' %(end_effect)
        outname_2='scaling_exp_idx_used_ee%i_power_CA.csv' % (end_effect)
        outname_3='scaling_exp_fit_ee%i_power_CA.csv' % (end_effect)

    np.savetxt('%s/%s' % (outdir, outname_1), c[0:8], delimiter=', ')
    np.savetxt('%s/%s' % (outdir, outname_2), c[8].transpose(), delimiter=', ')
    np.savetxt('%s/%s' % (outdir, outname_3), c[9].transpose(), delimiter=', ')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def run_motif_RG(CP, outdir, R1_idx, R2_idx):
    status_message('Motif RG [%i to %i]' %(R1_idx, R2_idx), outdir)
    MOTIF_RG = CP.get_radius_of_gyration(R1=R1_idx, R2=R2_idx)
    MEAN_MOTIF_RG = arrayfy(np.mean(MOTIF_RG))
    STD_MOTIF_RG  = arrayfy(np.std(MOTIF_RG))
    np.savetxt('%s/motif_%i_%i_RG.csv' % (outdir, R1_idx, R2_idx), MOTIF_RG, delimiter=', ')
    np.savetxt('%s/motif_%i_%i_RG_mean.csv' % (outdir, R1_idx, R2_idx), MEAN_MOTIF_RG, delimiter=', ')
    np.savetxt('%s/motif_%i_%i_RG_std.csv' % (outdir, R1_idx, R2_idx), STD_MOTIF_RG, delimiter=', ')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def run_SASA(CP, outdir, stride2use, probe_radius=0.14):
    status_message('SASA', outdir)

    # read in max values for sidechains and backbone, and construct a dictionary that allows
    # easy lookup for each residue type
    #max_sasa_vals = np.loadtxt('/work/alex/tools/ANALYZER/data/SASA_SUMMARY.csv',delimiter=',')
    max_sasa_vals = MAX_SASA_DATA
    AA_max={}
    for idx in range(0,20):
        AA = AALIST[idx]
        AA_max[AA] = list(max_sasa_vals[idx])

    # get AA sequence (note this includes caps as  '<' and '>' but SASA also gives you
    AA_seq = CP.get_amino_acid_sequence(oneletter=True)

    NCAP=False
    CCAP=False
    if AA_seq[0] == '<':
        NCAP=True

    if AA_seq[-1] == '>':
        CCAP=True


    (SASA, SASA_SC, SASA_BB)  = CP.get_all_SASA(stride=stride2use, mode='all',probe_radius=probe_radius)
    #SASA_BB = CP.get_all_SASA(stride=stride, mode='backbone')
    #SASA_SC = CP.get_all_SASA(stride=stride, mode='sidechain')

    if NCAP:
        SASA = SASA.transpose()[1:].transpose()
        AA_seq = AA_seq[1:]
    if CCAP:
        SASA = SASA.transpose()[:-1].transpose()
        AA_seq = AA_seq[:-1]

    # quick sanity check
    if SASA.shape[0] != SASA_BB.shape[0] or SASA.shape[0] != SASA_SC.shape[0]:
        raise AnalyzerException('Error - mismatch in SASA analysis... -> suggests a bug')

    if SASA.shape[1] != SASA_BB.shape[1] or SASA.shape[1] != SASA_SC.shape[1]:
        raise AnalyzerException('Error - mismatch in SASA analysis... -> suggests a bug')


    # build max_vectors - vectors that are length of sequence and have max possible SAASA
    # associated with each residue position
    BB_max=[]
    SC_max=[]
    for i in range(0,SASA.shape[1]):
        SC_max.append(AA_max[AA_seq[i]][0])
        BB_max.append(AA_max[AA_seq[i]][1])


    norm_bb = SASA_BB/BB_max
    norm_sc = SASA_SC/SC_max
    norm_sc[np.isnan(norm_sc)]=-1

    # absolute vals
    MEAN_SASA = np.mean(SASA,0)
    STD_SASA =  np.std(SASA,0)

    MEAN_SASA_BB = np.mean(SASA_BB,0)
    STD_SASA_BB =  np.std(SASA_BB,0)

    MEAN_SASA_SC = np.mean(SASA_SC,0)
    STD_SASA_SC =  np.std(SASA_SC,0)

    # norm vals
    MEAN_SASA_BB_norm = np.nanmean(norm_bb,0)
    STD_SASA_BB_norm =  np.nanstd(norm_bb,0)

    MEAN_SASA_SC_norm = np.nanmean(norm_sc,0)
    STD_SASA_SC_norm =  np.nanstd(norm_sc,0)

    # get amino acid sequence and

    if probe_radius == 0.14:
        np.savetxt('%s/SASA_mean.csv' % (outdir), MEAN_SASA, delimiter=', ')
        np.savetxt('%s/SASA_std.csv' % (outdir), STD_SASA, delimiter=', ')

        np.savetxt('%s/SASA_BB_mean.csv' % (outdir), MEAN_SASA_BB, delimiter=', ')
        np.savetxt('%s/SASA_BB_std.csv' % (outdir), STD_SASA_BB, delimiter=', ')

        np.savetxt('%s/SASA_SC_mean.csv' % (outdir), MEAN_SASA_SC, delimiter=', ')
        np.savetxt('%s/SASA_SC_std.csv' % (outdir), STD_SASA_SC, delimiter=', ')

        np.savetxt('%s/SASA_BB_mean_norm.csv' % (outdir), MEAN_SASA_BB_norm, delimiter=', ')
        np.savetxt('%s/SASA_BB_std_norm.csv' % (outdir), STD_SASA_BB_norm, delimiter=', ')

        np.savetxt('%s/SASA_SC_mean_norm.csv' % (outdir), MEAN_SASA_SC_norm, delimiter=', ')
        np.savetxt('%s/SASA_SC_std_norm.csv' % (outdir), STD_SASA_SC_norm, delimiter=', ')

    else:
        np.savetxt('%s/SASA_mean_radius_%2.2f.csv' % (outdir, probe_radius), MEAN_SASA, delimiter=', ')
        np.savetxt('%s/SASA_std_radius_%2.2f.csv' % (outdir, probe_radius), STD_SASA, delimiter=', ')

        np.savetxt('%s/SASA_BB_mean_radius_%2.2f.csv' % (outdir, probe_radius), MEAN_SASA_BB, delimiter=', ')
        np.savetxt('%s/SASA_BB_std_radius_%2.2f.csv' % (outdir, probe_radius), STD_SASA_BB, delimiter=', ')

        np.savetxt('%s/SASA_SC_mean_radius_%2.2f.csv' % (outdir, probe_radius), MEAN_SASA_SC, delimiter=', ')
        np.savetxt('%s/SASA_SC_std_radius_%2.2f.csv' % (outdir, probe_radius), STD_SASA_SC, delimiter=', ')


def run_DSSP_analysis_OLD(CP, outdir):

    dssp_data = md.compute_dssp(CP.traj)
    C_vector = []
    E_vector = []
    H_vector = []

    n_residues = CP.n_residues
    n_frames   = CP.n_frames

    for i in range(1,n_residues-1):
        C_vector.append(float(sum(dssp_data.transpose()[i] == 'C'))/n_frames)
        E_vector.append(float(sum(dssp_data.transpose()[i] == 'E'))/n_frames)
        H_vector.append(float(sum(dssp_data.transpose()[i] == 'H'))/n_frames)

    np.savetxt('%s/DSSP_H.csv' % (outdir), np.array(H_vector), delimiter=', ')
    np.savetxt('%s/DSSP_E.csv' % (outdir), np.array(E_vector), delimiter=', ')
    np.savetxt('%s/DSSP_C.csv' % (outdir), np.array(C_vector), delimiter=', ')


def run_DSSP_analysis(CP, outdir):

    dssp_out = CP.get_secondary_structure_DSSP()

    np.savetxt('%s/DSSP_H.csv' % (outdir), dssp_out[1], delimiter=', ')
    np.savetxt('%s/DSSP_E.csv' % (outdir), dssp_out[2], delimiter=', ')
    np.savetxt('%s/DSSP_C.csv' % (outdir), dssp_out[3], delimiter=', ')


def run_BBSEG_analysis(CP, outdir):
    bbseg_out = CP.get_secondary_structure_BBSEG()

    for i in range(0, 9):
        np.savetxt('%s/BBSEG_%i.csv' % (outdir,i), np.array(bbseg_out[i]), delimiter=', ')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def run_linear_heterogeneity(CP, outdir):
    LH = CP.get_local_heterogeneity(stride=20)
    np.savetxt('%s/linear_heterogeneity_mean.csv' % outdir, np.hstack((np.zeros(5), LH[0], np.zeros(5))), delimiter=', ')
    np.savetxt('%s/linear_heterogeneity_std.csv' % outdir, np.hstack((np.zeros(5), LH[1], np.zeros(5))), delimiter=', ')


def run_heterogeneity_analysis(CP, strideval, outdir):
    ALL_D = CP.get_D_vector(stride=strideval)
    mean_D  = arrayfy(np.mean(ALL_D))
    std_D   = arrayfy(np.std(ALL_D))
    np.savetxt('%s/D_vector.csv' % (outdir), ALL_D, delimiter=', ')
    np.savetxt('%s/D_mean.csv' % (outdir), mean_D, delimiter=', ')
    np.savetxt('%s/D_std.csv' % (outdir), std_D, delimiter=', ')


def run_cluster_analysis(CP, strideval, outdir):

    glob = CP.get_clusters(stride=strideval, n_clusters=10)
    np.savetxt('%s/cluster_size.csv' % (outdir), np.array(glob[0]), delimiter=', ')

    centroids = glob[3]

    c = centroids[0]
    t = glob[1][0]

    centroid_traj = t[c]
    for idx in range(1,10):
        c = centroids[idx]
        t = glob[1][idx]

        centroid_traj = centroid_traj + t[c]


    centroid_traj[0].save_pdb('cluster_centroid_traj.pdb')
    centroid_traj.save_xtc('cluster_centroid_traj.xtc')


def run_dihedral_extraction(CP, outdir):
    MEGA_PHI=[]
    MEGA_PSI=[]
    MEGA_OMEGA=[]
    for residue_index in range(1, CP.n_residues - 3):
        print(residue_index)
        PHI   = np.degrees(np.transpose(md.compute_phi(CP.traj)[1])[residue_index])
        PSI   = np.degrees(np.transpose(md.compute_psi(CP.traj)[1])[residue_index])
        OMEGA = np.degrees(np.transpose(md.compute_omega(CP.traj)[1])[residue_index])

        MEGA_PHI.append(PHI)
        MEGA_PSI.append(PSI)
        MEGA_OMEGA.append(OMEGA)


    np.savetxt('%s/PSI_matrix.csv' % (outdir), np.array(MEGA_PSI), delimiter=', ')
    np.savetxt('%s/PHI_matrix.csv' % (outdir), np.array(MEGA_PHI), delimiter=', ')
    np.savetxt('%s/OMEGA_matrix.csv' % (outdir), np.array(MEGA_OMEGA), delimiter=', ')


def run_angle_mutual_information(CP, outdir, angle_name):
    MIMatrix = CP.get_dihedral_mutual_information(angle_name=angle_name)
    np.savetxt('%s/%s_mutual_information.csv' % (outdir, angle_name), MIMatrix, delimiter=', ')
