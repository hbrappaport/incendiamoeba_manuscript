import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

plt.style.use('default')
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
df["Positive Surface Charge Fraction"] = df["num_exposed_pos"] / (df["num_exposed_pos"] + df["num_exposed_neg"])
df["Negative Surface Charge Fraction"] = df["num_exposed_neg"] / (df["num_exposed_pos"] + df["num_exposed_neg"])

df["Net Surface Charge"] = df["surface_net_charge"]

# Define charge parameters to analyze

name = "charge_analysis_main"
params = [
    "Net Surface Charge",
    "Positive Surface Charge Fraction",
    "Negative Surface Charge Fraction",
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
    count=0
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
                line_color = 'blue'
            elif change < 0:
                line_color = 'red'
            
            ax.plot([org_x_positions[organisms[0]], org_x_positions[organisms[1]]],
                   [val1, val2], color=line_color, alpha=0.6, linewidth=1.5, zorder=0)
            count +=1
    print(count)
            


    
    # Add significance indication based on FDR-corrected p-value
    if len(all_org1_values) > 0 and len(all_org2_values) > 0:

        
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
    ax.set_ylabel(param)
    ax.set_xlim(-0.25, 1.25)





# Remove extra subplots
for p in range(n_params, len(axes)):
    fig.delaxes(axes[p])

plt.tight_layout()
plt.savefig(f'{name}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{name}.svg',dpi = 300, bbox_inches = 'tight', pad_inches=0.1)


plt.close()
