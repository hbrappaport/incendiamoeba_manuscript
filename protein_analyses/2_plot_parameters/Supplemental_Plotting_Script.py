import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

"""
Author: Jess Niblo, jkniblo@syr.edu
Date: 2025-09-30

Purpose:
    Plot metrics from analyze_struc_params_single.py
"""
plt.style.use('default')

plt.rcParams['svg.fonttype'] = 'none' #set the value globally
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['xtick.minor.width'] = 2.5
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 2.5
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.minor.size'] = 8
plt.rcParams['ytick.minor.width'] = 2.5
plt.rcParams['font.size'] = 12
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['legend.title_fontsize'] = 12

exp_file = "Incendi_Vverm_proteinpredictions.csv"
df = pd.read_csv(exp_file)

df["residue_per_Rg"] = df["num_residues"] / df["Rg nm"]
df["sasa_per_residue"] = df["total_sasa"] / df["num_residues"]
df["disulfide_per_Rg"] = df["num_disulfides"] / df["Rg nm"]
df["salt_bridges_per_Rg"] = df["num_salt_bridge"] / df["Rg nm"]
df["hbond_per_Rg"] = df["num_hbonds"] / df["Rg nm"] 
df["pipi_per_Rg"] = df["num_pipi"] / df["Rg nm"] 
df["helix_per_Rg"] = df["helix_num"] / df["Rg nm"]
df["strand_per_Rg"] = df["strand_num"] / df["Rg nm"]
df["coil_per_Rg"] = df["coil_num"] / df["Rg nm"]
df["hbonds_per_helix"] = df["num_hbonds"] / df["helix_num"].replace(0, np.nan)
df["hbonds_per_strand"] = df["num_hbonds"] / df["strand_num"].replace(0, np.nan)
df["salt_bridges_per_helix"] = df["num_salt_bridge"] / df["helix_num"].replace(0, np.nan)
df["salt_bridges_per_strand"] = df["num_salt_bridge"] / df["strand_num"].replace(0, np.nan)
df["disulfide_per_SS"] = df["num_disulfides"] / (df["helix_num"] + df["strand_num"] + df["coil_num"])
df["frac_pos_exposed"] = df["num_exposed_pos"] / df["Total_pos_charge"]
df["frac_neg_exposed"] = df["num_exposed_neg"] / df["Total_neg_charge"]
df["surface_charge_per_SASA"] = (df["num_exposed_pos"] + df["num_exposed_neg"]) / df["total_sasa"]
df["frac_sasa_aromatic"] = (df["frac_sasa_TRP"] + df["frac_sasa_TYR"] + df["frac_sasa_PHE"]) / df["total_sasa"]


ylabels = {
    "residue_per_Rg": "Residues per\nRadius of Gyration",
    "disulfide_per_Rg": "Disulfide Bonds per\nRadius of Gyration",
    "salt_bridges_per_Rg": "Salt Bridges per\nRadius of Gyration",
    "hbond_per_Rg": "Hydrogen Bonds per\nRadius of Gyration",
    "pipi_per_Rg": "π–π Interactions per\nRadius of Gyration",
    "helix_per_Rg": "Helices per\nRadius of Gyration",
    "strand_per_Rg": "Strands per\nRadius of Gyration",
    "coil_per_Rg": "Coils per\nRadius of Gyration",
    "hbonds_per_helix": "Hydrogen Bonds per Helix",
    "hbonds_per_strand": "Hydrogen Bonds per Strand",
    "salt_bridges_per_helix": "Salt Bridges per Helix",
    "salt_bridges_per_strand": "Salt Bridges per Strand",
    "disulfide_per_SS": "Disulfide Bonds per\nSecondary Structure",
    "Total_ratio_neg_pos": "Ratio of Negative to\nPositive Residues",
    "Theoretical_pI": "Theoretical Isoelectric Point",
    "total_sasa": "Total SASA",
    "surface_charge_per_SASA": "Total Surface Charge per SASA",
    "frac_pos_exposed": "Fraction of Exposed\nPositive Residues",
    "frac_neg_exposed": "Fraction of Exposed\nNegative Residues",
    "frac_sasa_hydrophobic": "Hydrophobic SASA Fraction",
    "frac_sasa_aromatic": "Aromatic SASA Fraction",
    "frac_sasa_polar": "Polar SASA Fraction",
    "frac_sasa_charged": "Charged SASA Fraction",
    "frac_sasa_pos_charge": "Positively Charged\nSASA Fraction",
    "frac_sasa_neg_charge": "Negatively Charged\nSASA Fraction"
}


name = "mega_plot_FDR_SI"
params = [

    #! structure 
    "residue_per_Rg",

    #! bonds per size 
    "disulfide_per_Rg",
    "salt_bridges_per_Rg",
    "hbond_per_Rg",
    "pipi_per_Rg",

    #! secondary structure and interactions per SS 
    "helix_per_Rg",
    "strand_per_Rg",
    "coil_per_Rg",
    "hbonds_per_helix",
    "hbonds_per_strand",
    "salt_bridges_per_helix",
    "salt_bridges_per_strand",
    "disulfide_per_SS",

    #! charge distributions - Global
    "Total_ratio_neg_pos",
    "Theoretical_pI",

    #! charge distributions - Surface
    "total_sasa",
    "surface_charge_per_SASA",
    "frac_pos_exposed",
    "frac_neg_exposed",
   
    #! SASA by amino acid group
    "frac_sasa_hydrophobic",
    "frac_sasa_aromatic",
    "frac_sasa_polar",
    "frac_sasa_charged",
    "frac_sasa_pos_charge",
    "frac_sasa_neg_charge",

]


unique_proteins = df["Protein"].unique()
organisms = ["Incendiamoeba", "Vermamoeba"]
organisms_label = ["I. cascadae", "V. vermiformis"]


n_params = len(params)
n_cols = 5
n_rows = int(np.ceil(n_params / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(3*n_cols, 3*n_rows))
fig.subplots_adjust(wspace=0.05, hspace=0.25)
axes = axes.flatten()

p_values = []
param_names = []
test_results = []

# First pass: collect all p-values
for p, param in enumerate(params):
    
    all_org1_values = []
    all_org2_values = []
    
    for protein in unique_proteins:
        org_values = {}
        
        for i, org in enumerate(organisms):
            subset = df[(df["Protein"] == protein) & (df["Organism"] == org)]
            value = subset[param].values
            
            if len(value) == 0 or np.all(np.isnan(value)):
                mean_value = np.nan
            else:
                mean_value = np.nanmean(value)
            
            org_values[org] = mean_value
        
        val1 = org_values[organisms[0]]
        val2 = org_values[organisms[1]]
        
        if not (np.isnan(val1) or np.isnan(val2)):
            all_org1_values.append(val1)
            all_org2_values.append(val2)
    
    # Perform paired t-test
    if len(all_org1_values) > 0 and len(all_org2_values) > 0:
        t_stat, p_value = stats.ttest_rel(all_org1_values, all_org2_values)
        p_values.append(p_value)
        param_names.append(param)
        test_results.append({
            'org1_values': all_org1_values,
            'org2_values': all_org2_values,
            'p_value': p_value
        })
    else:
        p_values.append(1.0)
        param_names.append(param)
        test_results.append({
            'org1_values': [],
            'org2_values': [],
            'p_value': 1.0
        })

# Apply FDR correction using Benjamini-Hochberg procedure
fdr_corrected = false_discovery_control(p_values, method='bh')

# Create plots with corrected p-values
for p, param in enumerate(params):
    ax = axes[p]
    
    all_org1_values = test_results[p]['org1_values']
    all_org2_values = test_results[p]['org2_values']
    raw_p_value = test_results[p]['p_value']
    corrected_p_value = fdr_corrected[p]
    
    # Plot data points and connections
    for protein in unique_proteins:
 
        org_values = {}
        org_x_positions = {}
        
        for i, org in enumerate(organisms):
            subset = df[(df["Protein"] == protein) & (df["Organism"] == org)]
            value = subset[param].values
            
            if len(value) == 0 or np.all(np.isnan(value)):
                mean_value = np.nan
            else:
                mean_value = np.nanmean(value)
            
            org_values[org] = mean_value
            x_jitter = np.random.normal(0, 0.03)
            x_positions = i + x_jitter
            org_x_positions[org] = x_positions
            
            ax.scatter(x_positions, mean_value, color='white', edgecolors='black', 
                      linewidths=1, zorder=5, s=30)
        
        val1 = org_values[organisms[0]]
        val2 = org_values[organisms[1]]
        
        if not (np.isnan(val1) or np.isnan(val2)):
            change = val2 - val1
            if change > 0: 
                line_color = '#3A6AA1'
            elif change < 0:
                line_color = '#C1532D'
            
            ax.plot([org_x_positions[organisms[0]], org_x_positions[organisms[1]]],
                   [val1, val2], color=line_color, alpha=0.6, linewidth=1.5, zorder=0)

    
    # Add significance indication based on FDR-corrected p-value
    if len(all_org1_values) > 0 and len(all_org2_values) > 0:
        if raw_p_value < 0.0005:
            sig_label = "***"
        elif raw_p_value < 0.005:
            sig_label = "**"
        elif raw_p_value < 0.05:
            sig_label = "*"
        else:
            sig_label = "ns"
        
        # Color background by FDR-corrected p-value
        if corrected_p_value >= 0.05:
            bg_color = 'white'
        else:
            bg_color = 'lightgray'
        
        ax.set_facecolor(bg_color)
        
        # Add significance bracket and label
        y_max = max(max(all_org1_values), max(all_org2_values))
        y_min = min(min(all_org1_values), min(all_org2_values))

        y_range = y_max - y_min
        if y_range == 0 or np.isclose(y_range, 0):
            y_range = 1.0  # fallback so scaling is reasonable

        bracket_height = y_max + 0.05 * y_range
        ax.set_ylim(y_min - 0.1 * y_range, bracket_height + 0.2 * y_range)



        ax.plot([0, 0, 1, 1], 
               [bracket_height, bracket_height + 0.02 * y_range, 
                bracket_height + 0.02 * y_range, bracket_height], 
               'k-', linewidth=1.5)
        

        ax.text(0.5, bracket_height , 
            f"FDR adj. P={corrected_p_value:.1}", 
            ha='center', va='bottom', fontsize=11)
        

        
        
    else:
        ax.set_facecolor('lightgray')
    
    # Format subplot
    ax.set_xticks([0, 1])
    ax.set_xticklabels(organisms_label, fontstyle='italic')
    ax.set_ylabel(ylabels[param])

    ax.set_xlim(-0.25, 1.25)





# Remove extra subplots
for p in range(n_params, len(axes)):
    fig.delaxes(axes[p])

plt.tight_layout()
plt.savefig(f'{name}_3.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{name}_3.svg',dpi = 300, bbox_inches = 'tight', pad_inches=0.1)


plt.close()
