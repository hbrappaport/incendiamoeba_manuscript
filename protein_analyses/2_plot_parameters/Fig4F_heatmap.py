import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap

plt.style.use('default')
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 8

SPECIES = [
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Dictyostelium_discoideum_comparison.csv",  "D. discoideum",    27,  "Meso. Amoeba"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/vvermi_comparison.csv",                    "V. vermiformis",   45,  "Meso./Thermo Tol. Amoeba"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Porphyra_umbilicalis_comparison.csv",      "P. umbilicalis",   34,  "Meso. Algea"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Cyanidioschyzon_merolae_comparison.csv",   "C. merolae",       56,  "Thermo. Algae"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Talaromyces_amestolkiae_comparison.csv",   "T. amestolkiae",   37,  "Meso. Fungi"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Aspergillus_sergii_comparison.csv",        "A. sergii",        40,  "Meso. Fungi"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Chaetomium_thermophilum_comparison.csv",   "C. thermophilum",  60,  "Thermo. Fungi"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Thermomyces_lanuginosus_comparison.csv",   "T. lanuginosus",   63,  "Thermo. Fungi"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Herpetosiphon_gulosus_comparison.csv",     "H. gulosus",       37,  "Meso. Bacteria (Chloroflexi)"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Caldilinea_aerophila_comparison.csv",      "C. aerophila",     65,  "Thermo. Bacteria (Chloroflexi)"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Deinococcus_ficus_comparison.csv",         "D. ficus",         35,  "Meso. Bacteria (Deinococcota)"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Thermus_aquaticus_comparison.csv",         "T. aquaticus",     80,  "Thermo. Bacteria (Deinococcota)"),
    ("../1_analyze_pdb_structures/structural_parameters_raw_data/Pyrobaculum_ferrireducens_comparison.csv", "P. ferrireducens", 100, "Thermo. Archaea"),
]

DOMAIN_SPLIT_AFTER = 7
OUT_NAME           = "Fig4F_heatmap"

INCENDI_DISORDER_CSV = "/home/jkniblo/Ameoba/quantify_disorder/incendi_500_metapredictv3_Mar2026.csv"

INCENDI_FASTA = "fastas/incendi_500.fasta"

FASTA_FILES = {
    "D. discoideum":    "fastas/Dictyostelium_discoideum.fasta",
    "V. vermiformis":   "fastas/vvermi.fasta",
    "P. umbilicalis":   "fastas/Porphyra_umbilicalis.fasta",
    "C. merolae":       "fastas/Cyanidioschyzon_merolae.fasta",
    "T. amestolkiae":   "fastas/Talaromyces_amestolkiae.fasta",
    "A. sergii":        "fastas/Aspergillus_sergii.fasta",
    "C. thermophilum":  "fastas/Chaetomium_thermophilum.fasta",
    "T. lanuginosus":   "fastas/Thermomyces_lanuginosus.fasta",
    "H. gulosus":       "fastas/Herpetosiphon_gulosus.fasta",
    "C. aerophila":     "fastas/Caldilinea_aerophila.fasta",
    "D. ficus":         "fastas/Deinococcus_ficus.fasta",
    "T. aquaticus":     "fastas/Thermus_aquaticus.fasta",
    "P. ferrireducens": "fastas/Pyrobaculum_ferrireducens.fasta",
}

nonpolar_uncharged = ["ALA", "VAL", "ILE", "LEU", "MET"]

# ─── DERIVED METRICS ──────────────────────────────────────────────────────────
def add_derived_metrics(df):
    for prefix in ["incendi_", "ortholog_"]:
        p = prefix
        npu_cols = [f"{p}frac_sasa_{aa}" for aa in nonpolar_uncharged if f"{p}frac_sasa_{aa}" in df.columns]
        total_surface = df[f"{p}num_exposed_pos"] + df[f"{p}num_exposed_neg"]

        df[f"{p}salt_bridges_per_Rg"] = df[f"{p}num_salt_bridges"] / df[f"{p}rg_nm"]
        df[f"{p}frac_sasa_nonpolar_uncharged"] = df[npu_cols].sum(axis=1)
        df[f"{p}Positive_Surface_Charge_Fraction"] = df[f"{p}num_exposed_pos"] / total_surface.replace(0, np.nan)
        df[f"{p}Negative_Surface_Charge_Fraction"] = df[f"{p}num_exposed_neg"] / total_surface.replace(0, np.nan)

    return df

def extract_protein_id(full_name):
    parts = str(full_name).split('|')
    if len(parts) >= 2:
        return parts[1].strip()
    return str(full_name).split()[0].strip()

# ─── PARAMS ───────────────────────────────────────────────────────────────────
params = [
    "frac_sasa_SER",
    "frac_sasa_hydrophobic",
    "frac_sasa_nonpolar_uncharged",
    "Positive_Surface_Charge_Fraction", "Negative_Surface_Charge_Fraction",
    "frac_sasa_LYS", "frac_sasa_ARG",
    "salt_bridges_per_Rg",

]

param_labels = {
    "frac_sasa_hydrophobic":            "Hydrophobic SASA Frac.",
    "frac_sasa_SER":                    "Serine SASA Frac.",
    "frac_sasa_LYS":                    "Lysine SASA Frac.",
    "frac_sasa_ARG":                    "Arginine SASA Frac.",
    "salt_bridges_per_Rg":              "Salt Bridges per $R_g$",
    "Positive_Surface_Charge_Fraction": '+ Surface Charge Frac.',
    "Negative_Surface_Charge_Fraction": '- Surface Charge Frac.',
    "frac_sasa_nonpolar_uncharged":     "Nonpolar SASA Frac.",
    
}


# Store paired arrays (incendi, ortholog) per param per species for paired t-test
species_pairs  = {param: {} for param in params}
species_labels = []
species_ns     = {}
species_temp   = {}
species_class  = {}

for csv_path, label, temp, classification in SPECIES:
    path = Path(csv_path)
    if not path.exists():
        print(f"  SKIP {label}: file not found ({csv_path})")
        continue

    df = pd.read_csv(path)
    df = add_derived_metrics(df)

    species_labels.append(label)
    species_ns[label]    = len(df)
    species_temp[label]  = temp
    species_class[label] = classification


    # store paired arrays 
    for param in params:
        col_i = f"incendi_{param}"
        col_o = f"ortholog_{param}"
        if col_i not in df.columns or col_o not in df.columns:
            continue
        iv   = df[col_i].values.astype(float)
        ov   = df[col_o].values.astype(float)
        mask = ~(np.isnan(iv) | np.isnan(ov))
        species_pairs[param][label] = (iv[mask], ov[mask])


# Drop params with no data in any species
params = [
    p for p in params
    if any(len(v[0]) > 0 for v in species_pairs[p].values())
]
print(f"\n{len(params)} parameters, {len(species_labels)} species")

#Pair t-test FDR corrected 
all_pvals  = []
pval_index = []
for param in params:
    for label in species_labels:
        pair = species_pairs[param].get(label, (np.array([]), np.array([])))
        iv_p, ov_p = pair
        if len(iv_p) >= 2:
            _, p = stats.ttest_rel(iv_p, ov_p)
            p = 1.0 if (np.isnan(p) or np.isinf(p)) else p
        else:
            p = 1.0
        all_pvals.append(p)
        pval_index.append((param, label))

fdr_corrected = false_discovery_control(all_pvals, method='bh')
fdr_map = {key: fdr for key, fdr in zip(pval_index, fdr_corrected)}

n_params  = len(params)
n_species = len(species_labels)

median_matrix = np.full((n_species, n_params), np.nan)
sig_matrix    = np.zeros((n_species, n_params), dtype=bool)

for j, param in enumerate(params):
    for i, label in enumerate(species_labels):
        pair = species_pairs[param].get(label, (np.array([]), np.array([])))
        delta = pair[0] - pair[1]
        if len(delta) > 0:
            median_matrix[i, j] = np.median(delta)
        sig_matrix[i, j] = fdr_map.get((param, label), 1.0) < 0.01

# normalization within each column 
norm_matrix = np.full_like(median_matrix, np.nan)
for j in range(n_params):
    col = median_matrix[:, j]
    col_abs_max = np.nanmax(np.abs(col))
    if col_abs_max > 0:
        norm_matrix[:, j] = col / col_abs_max
    else:
        norm_matrix[:, j] = 0.0



cell_w, cell_h = 0.55, 0.55
fig_w = max(8, n_params  * cell_w + 4)
fig_h = max(2, n_species * cell_h + 2)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

cmap = LinearSegmentedColormap.from_list(
    'custom',
    ['#7b5ea7', '#fafaf8', '#d4a843'],
    N=256
)
im = ax.imshow(
    norm_matrix,
    cmap=cmap,
    aspect='auto',
    vmin=-1, vmax=1,
    interpolation='nearest',
)

#Significance Stars 
for i in range(n_species):
    for j in range(n_params):
        if sig_matrix[i, j]:
            ax.text(j, i, '★', ha='center', va='center',
                    fontsize=16, color='black', zorder=3)

x_labels = [param_labels.get(p, p) for p in params]
ax.set_xticks(range(n_params))
ax.set_xticklabels(x_labels, fontsize=16, rotation=90)

ax.set_yticks(range(n_species))
y_labels = []
for label in species_labels:
    italic_label   = label.replace(' ', '\\ ')
    classification = species_class[label]
    y_labels.append(f"$\\it{{{italic_label}}}$") 
ax.set_yticklabels(y_labels, fontsize=16)
for ytick, label in zip(ax.get_yticklabels(), species_labels):
    ytick.set_color('k')

ax.set_xticks(np.arange(n_params) - 0.5, minor=True)
ax.set_yticks(np.arange(n_species) - 0.5, minor=True)
ax.grid(which='minor', color='white', linewidth=0.8)
ax.tick_params(which='minor', length=0)
ax.tick_params(which='major', length=6, width=1.5)

##Domain divider
if DOMAIN_SPLIT_AFTER is not None:
    ax.axhline(
        y=DOMAIN_SPLIT_AFTER + 0.5,
        color='black',
        linewidth=2.0,
        linestyle='--',
        zorder=4,
    )

cbar = fig.colorbar(im, ax=ax, pad=0.01)
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels(
    ['Higher in\nOrtholog', 'No\ndifference', 'Higher in\nIncendi'],
    fontsize=14,
)
cbar.ax.tick_params(length=6, width=1.5)

fig.tight_layout()
plt.savefig(f"{OUT_NAME}.pdf", dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.close()
