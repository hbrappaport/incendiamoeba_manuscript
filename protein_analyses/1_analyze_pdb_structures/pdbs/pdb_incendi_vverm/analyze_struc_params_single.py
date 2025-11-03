import sys
import numpy as np 
import mdtraj as md
import MDAnalysis as mda
import os 
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import pandas as pd
from MDAnalysis.analysis import contacts

"""
Author: Jess Niblo, jkniblo@syr.edu
Date: 2025-08-27

Purpose:
    Analysis of PDB structures for strucrual properities, surface properties, and molecular interactions.

Dependencies:
    see StrucParams.yml (run conda env create -f StrucParams.yml)

How to use: 
    1. place all pdb files of interest in the same directory as this script
    2. fill in list of pdbs to be analyzed on line 248
    3. Verify cutoffs are appropriate for your system -- lines 54, 94, 179, 203, 239-244
    4. the results will be saved in the output directory 
"""


hydrophobic = ['ALA','VAL','ILE','LEU','MET','PHE','PRO','TRP']
polar       = ['SER','THR','ASN','GLN','TYR','CYS','GLY']
charged     = ['ASP','GLU','LYS','ARG']
positive = ['ARG', 'LYS'] 
negative = ['ASP', 'GLU']


def contacts_within_cutoff(u, group_a, group_b, cutoff):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= cutoff
        n_contacts = contacts.contact_matrix(dist, cutoff).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)


def find_salt_bridges(struc_mdA): 
    ## MDAnalysis
    ##following prodecure here: https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html
    sel_basic = "(resname ARG LYS) and (name NH* NZ)"
    sel_acidic = "(resname ASP GLU) and (name OE* OD*)"
    acidic = struc_mdA.select_atoms(sel_acidic)
    basic = struc_mdA.select_atoms(sel_basic)
    ca = contacts_within_cutoff(struc_mdA, acidic, basic, 4.5)
    num_bridges = ca[0][1] ##first column = frame, second column = number of contacts

    return num_bridges

def find_disulfides(pdb_file, cutoff):
    # select all sulfur atoms from cysteines
    struc_md = md.load_pdb(pdb_file)
    sulfurs = struc_md.topology.select('resname CYS and name SG')
    positions = struc_md.xyz[0, sulfurs, :]  # single frame

    n = len(sulfurs)
    disulfides = 0
    for i in range(n):
        for j in range(i+1, n):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist <= cutoff:  # Å
                disulfides += 1 #.append((sulfurs[i], sulfurs[j]))
    return disulfides


def analyze_hbonds(PDB_file): 
    h_filename = add_hydrogens(PDB_file)
    mdt_hy = md.load_pdb(f"{h_filename}.pdb")
    hbonds = md.baker_hubbard(mdt_hy,periodic=False) #hbonds 

    num_hbonds = len(hbonds) 
    
    return num_hbonds

def find_pi_pi_interactions(struc_mdA):
    ## MDAnalysis 
    aromatic_residues = struc_mdA.select_atoms("resname PHE or resname TYR or resname TRP").residues
    ring_groups = [res.atoms.select_atoms("not name H*") for res in aromatic_residues]



    counts = 0
    for i in range(len(ring_groups)):
        for j in range(i+1, len(ring_groups)):
            pair_contacts = contacts_within_cutoff(struc_mdA, ring_groups[i], ring_groups[j], 4.5)

            counts += pair_contacts[0][1]  ##first column = frame, second column = number of contacts
    return counts


def add_hydrogens(PDB_file): 
    base = os.path.splitext(os.path.basename(PDB_file))[0]


    fixer = PDBFixer(filename=PDB_file)
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(f"{base}_withH.pdb",'w'))

    return f"{base}_withH"

def analyze_geometry(struc): 

    num_chains = struc.topology.n_chains
    num_residues = struc.topology.n_residues
    rg = md.compute_rg(struc)

    ##secondary structure content 
    dssp = md.compute_dssp(struc, simplified=True)
    codes = dssp[0]  # single frame = one row

    #! uncomment comment to get a fraction
    helix_num = np.sum(codes == 'H') #/ len(codes) 
    strand_num = np.sum(codes == 'E') #/ len(codes)
    coil_num  = np.sum(codes == 'C') #/ len(codes)
    
    return num_chains, num_residues, rg[0], helix_num, strand_num, coil_num

def analyze_surface(struc_md): 
    non_h_atoms = [atom.index for atom in struc_md.topology.atoms if atom.element.symbol != 'H']
    traj_no_h = struc_md.atom_slice(non_h_atoms)

    sasa = md.shrake_rupley(traj_no_h, probe_radius=1.4) ##this on it's own is the SASA for each AA without hydrogen
    sasa *= 100 #convert from nm^2 to A^2
    total_sasa = sasa.sum(axis=1)

    ##get AA name for each atom 
    residue_name_per_atom = np.array([a.residue.name for a in traj_no_h.topology.atoms]) 
    residue_names = np.array([res.name for res in traj_no_h.topology.residues])


    ##Determine if each AA is hydrophobic, polar, or charged
    is_hydro = np.array([ (name in hydrophobic) for name in residue_name_per_atom ])
    is_polar = np.array([ (name in polar) for name in residue_name_per_atom ])
    is_charged = np.array([ (name in charged) for name in residue_name_per_atom ])
    is_positive = np.array([ (name in positive) for name in residue_name_per_atom ])
    is_negative = np.array([ (name in negative) for name in residue_name_per_atom ])

    ##Get SASA for each group
    hydro_val = sasa[:, is_hydro].sum(axis=1)
    polar_val = sasa[:, is_polar].sum(axis=1)
    charged_val = sasa[:, is_charged].sum(axis=1)
    pos_sasa = sasa[:, is_positive].sum(axis=1)
    neg_sasa = sasa[:, is_negative].sum(axis=1)


    ## Group SASA/Total SASA = fraction contribution of each group 
    hydro_frac = hydro_val / total_sasa
    polar_frac = polar_val / total_sasa
    charged_frac = charged_val / total_sasa
    frac_pos = pos_sasa / total_sasa
    frac_neg = neg_sasa / total_sasa

     # Fraction of SASA per residue type
    frac_sasa_per_residue = {}
    for res_type in set(residue_names):
        # select atoms belonging to residues of this type
        atom_idxs = [i for i, a in enumerate(residue_name_per_atom) if a == res_type]
        frac_sasa_per_residue[res_type] = sasa[:, atom_idxs].sum() / total_sasa

    #TODO : hydrophobic/charged patch analysis? 

    return  total_sasa[0], hydro_frac[0], polar_frac[0], charged_frac[0], frac_pos[0], frac_neg[0], frac_sasa_per_residue


def analyze_interactions(PDB_file, struc_mdA):
     
    num_hbonds = analyze_hbonds(PDB_file)

    num_salt_bridge = find_salt_bridges(struc_mdA) #salt bridges
    num_disulfides = find_disulfides(PDB_file, 2.2)
    num_pipi = find_pi_pi_interactions(struc_mdA) #pipi bonds 

    return num_hbonds, num_salt_bridge, num_pipi, num_disulfides


def analyze_electrostats(struc_mdT):
    residue_names = np.array([res.name for res in struc_mdT.topology.residues])
    
    ## Net Charge 
    pos_charge = np.sum(residue_names == "LYS") + np.sum(residue_names == "ARG")
    neg_charge = np.sum(residue_names == "ASP") + np.sum(residue_names == "GLU")
    net_charge = pos_charge - neg_charge
    

    ## Surface charge 

    non_h_atoms = [atom.index for atom in struc_mdT.topology.atoms if atom.element.symbol != 'H']
    traj_no_h = struc_mdT.atom_slice(non_h_atoms)


    sasa = md.shrake_rupley(traj_no_h, probe_radius=1.4, mode='residue')
    sasa *= 100 #convert from nm^2 to A^2
    residue_names = np.array([res.name for res in traj_no_h.topology.residues])
    exposure_threshold = 20.0 # A^2
    exposed_residues = sasa[0] > exposure_threshold

    surface_positive = 0
    surface_negative = 0

    for i, (res_name, is_exposed) in enumerate(zip(residue_names, exposed_residues)):
        if is_exposed:
            if res_name in ['LYS', 'ARG']:
                surface_positive += 1
            elif res_name in ['ASP', 'GLU']:
                surface_negative += 1
    
    surface_net_charge = surface_positive - surface_negative



    #Degree of exposure 
    #https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635
    max_sasa = {
        'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167,
        'GLU': 223, 'GLN': 225, 'GLY': 104, 'HIS': 224, 'ILE': 197,
        'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240, 'PRO': 159,
        'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174
    }

    buried = 0 
    intermediate = 0 
    exposed = 0 

    rsa_values = []
    for i, res_name in enumerate(residue_names):
        if res_name in max_sasa:
            rsa = sasa[0][i] / max_sasa[res_name]
            rsa_values.append(rsa)

            if rsa < 0.2:
                buried += 1
            elif 0.2 <= rsa < 0.5:
                intermediate += 1
            elif rsa >= 0.5:
                exposed += 1

    highly_exposed_frac = exposed / len(rsa_values)


    



    #! uncomment if avg values are needed
    # mean_exposure = np.mean(rsa_values)
    # median_exposure = np.median(rsa_values)
    # std_exposure = np.std(rsa_values)


    return net_charge, surface_positive, surface_negative, surface_net_charge, highly_exposed_frac

#TODO: Add in any pdb file names *without* extension that you want generated within csv file 
PDB_of_interest=['vverm_actin_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000', 'incendi_actin_461ef_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000',
'incendi_hsp20_bbabc_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000', 'vverm_hsp20_56e89_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000',
'incendi_steroldemethyl_8fbbf_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000', 'vverm_steroldemethyl_d91f2_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000',
'incendi_histoneh3_176c8_0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000', 'vverm_histoneh3_10146_0_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000',
'incendi_ubiquitinc_77f37_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000', 'incendi_dihydropteridine_10c11_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000',
'vverm_dihydropteridine_c3851_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000', 'vverm_ubiquitinc_33d34_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000',
'incendi_hsp90_12f90_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000', 'vverm_hsp90_8035d_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000',
'incendi_phosphomethyl_cb6f3_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000', 'incendi_rhofamily_9cfa5_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000',
'vverm_phosphomethyl_d8d9b_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000', 'vverm_rhofamily_f4c3b_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000',
'incendi_synaptotagmin_33c10_0_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000', 'vverm_synaptotagmin_49cc8_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000',
'incendi_hprt1_233cd_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000', 'vverm_hprt1_53237_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000',
'incendi_L18Ae_25e91_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000', 'vverm_L18Ae_e757e_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000',
'incendi_coilhelix_9a2df_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000', 'vverm_coilhelix_4759f_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000',
'incendi_cathepsin_5bf13_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000', 'vverm_cathepsin_8df99_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000',
'incendi_spo11_59e16_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000', 'vverm_spo11_2c27b_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000',
'incendi_ubiquitinlike_92337_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000', 'vverm_ubiquitinlike_6a013_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000',
'incendi_coactosin_77bbb_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000', 'vverm_coactosin_f9894_unrelaxed_rank_004_alphafold2_ptm_model_1_seed_000',
'incendi_histoneh2b_e6657_0_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000', 'vverm_histoneh2b_842d9_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000',
'incendi_gpat3_d2019_unrelaxed_rank_002_alphafold2_ptm_model_1_seed_000', 'vverm_gpat3_fd676_0_unrelaxed_rank_003_alphafold2_ptm_model_1_seed_000',
'incendi_lysocarboxy_435d1_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000', 'vverm_lysocarboxy_e9f13_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000']

results = []

for name in PDB_of_interest: 

    PDB_file = f"{name}.pdb"
    print(f"Analyzing {PDB_file}")

    #load pdb files into mdTraj or mdAnalysis
    mdA_u = mda.Universe(f"{PDB_file}")
    mdT_u = md.load_pdb(f"{PDB_file}")

    #! assumed that HIS was neutral 
    num_chains, num_residues, rg, helix_num, strand_num, coil_num = analyze_geometry(mdT_u)
    total_sasa, hydro_frac, polar_frac, charged_frac, frac_pos, frac_neg, frac_sasa_per_residue = analyze_surface(mdT_u)
    num_hbonds, num_salt_bridge, num_pipi, num_disulfides = analyze_interactions(PDB_file, mdA_u)
    net_charge, surface_positive, surface_negative, surface_net_charge, highly_exposed_frac = analyze_electrostats(mdT_u)


    result = {
        "Name": name,
        "num_chains": num_chains,
        "num_residues": num_residues,
        "Rg nm": rg,  # Rg is an array, take first element for single-frame PDB
        "helix_num": helix_num,
        "strand_num": strand_num,
        "coil_num": coil_num,
        "total_sasa A^2": total_sasa,  # sum across atoms for single frame
        "frac_sasa_hydrophobic": hydro_frac,
        "frac_sasa_polar": polar_frac,
        "frac_sasa_charged": charged_frac,
        "frac_sasa_pos_charge": frac_pos,
        "frac_sasa_neg_charge": frac_neg,
        "num_hbonds": num_hbonds,
        "num_salt_bridge": num_salt_bridge,
        "num_disulfides": num_disulfides,
        "num_pipi": num_pipi,
        "net_charge": net_charge,
        "num_exposed_pos": surface_positive,
        "num_exposed_neg": surface_negative,
        "surface_net_charge": surface_net_charge,
        "highly_exposed_frac": highly_exposed_frac
    }

    for res, frac in frac_sasa_per_residue.items():
        result[f"frac_sasa_{res}"] = frac[0]

    

    
    
    results.append(result)

df = pd.DataFrame(results)
df.to_csv("analysis_results.csv", index=False)

