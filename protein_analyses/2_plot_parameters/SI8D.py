import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

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
OUT_NAME           = "SI8D"

INCENDI_DISORDER_CSV = "/home/jkniblo/Ameoba/quantify_disorder/incendi_500_metapredictv3_Mar2026.csv"

DISORDER_FILES = {
    "D. discoideum":    "metapredict_output/Dictyostelium_discoideum_metapredictv3_Mar2026.csv",
    "V. vermiformis":   "metapredict_output/vvermi_metapredictv3_Mar2026.csv",
    "P. umbilicalis":   "metapredict_output/Porphyra_umbilicalis_metapredictv3_Mar2026.csv",
    "C. merolae":       "metapredict_output/Cyanidioschyzon_merolae_metapredictv3_Mar2026.csv",
    "T. amestolkiae":   "metapredict_output/Talaromyces_amestolkiae_metapredictv3_Mar2026.csv",
    "A. sergii":        "metapredict_output/Aspergillus_sergii_metapredictv3_Mar2026.csv",
    "C. thermophilum":  "metapredict_output/Chaetomium_thermophilum_metapredictv3_Mar2026.csv",
    "T. lanuginosus":   "metapredict_output/Thermomyces_lanuginosus_metapredictv3_Mar2026.csv",
    "H. gulosus":       "metapredict_output/Herpetosiphon_gulosus_metapredictv3_Mar2026.csv",
    "C. aerophila":     "metapredict_output/Caldilinea_aerophila_metapredictv3_Mar2026.csv",
    "D. ficus":         "metapredict_output/Deinococcus_ficus_metapredictv3_Mar2026.csv",
    "T. aquaticus":     "metapredict_output/Thermus_aquaticus_metapredictv3_Mar2026.csv",
    "P. ferrireducens": "metapredict_output/Pyrobaculum_ferrireducens_metapredictv3_Mar2026.csv",
}
DISORDER_THRESHOLD = 0.5

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
        df[f"{p}hbond_per_Rg"]        = df[f"{p}num_hbonds"]       / df[f"{p}rg_nm"]
        df[f"{p}salt_bridges_per_Rg"] = df[f"{p}num_salt_bridges"] / df[f"{p}rg_nm"]
        df[f"{p}pipi_per_Rg"]         = df[f"{p}num_pipi"]         / df[f"{p}rg_nm"]
        df[f"{p}disulfide_per_Rg"]    = df[f"{p}num_disulfides"]   / df[f"{p}rg_nm"]
        df[f"{p}helix_per_Rg"]        = df[f"{p}helix_num"]        / df[f"{p}rg_nm"]
        df[f"{p}strand_per_Rg"]       = df[f"{p}strand_num"]       / df[f"{p}rg_nm"]
        df[f"{p}coil_per_Rg"]         = df[f"{p}coil_num"]         / df[f"{p}rg_nm"]

        npu_cols = [f"{p}frac_sasa_{aa}" for aa in nonpolar_uncharged if f"{p}frac_sasa_{aa}" in df.columns]
        df[f"{p}frac_sasa_nonpolar_uncharged"] = df[npu_cols].sum(axis=1)
        df[f"{p}ionic_vs_hbond"] = df[f"{p}num_salt_bridges"] / df[f"{p}num_hbonds"].replace(0, np.nan)
        total_surface = df[f"{p}num_exposed_pos"] + df[f"{p}num_exposed_neg"]
        df[f"{p}Positive_Surface_Charge_Fraction"] = df[f"{p}num_exposed_pos"] / total_surface.replace(0, np.nan)
        df[f"{p}Negative_Surface_Charge_Fraction"] = df[f"{p}num_exposed_neg"] / total_surface.replace(0, np.nan)
        aro_cols = [f"{p}frac_sasa_TRP", f"{p}frac_sasa_TYR", f"{p}frac_sasa_PHE"]
        present  = [c for c in aro_cols if c in df.columns]
        if present:
            df[f"{p}frac_sasa_aromatic"] = df[present].sum(axis=1)
    return df

def parse_disorder_scores(score_str):
    try:
        cleaned = (str(score_str)
                   .replace('\n', ' ')
                   .replace('[', '')
                   .replace(']', '')
                   .strip())
        return np.array([float(x) for x in cleaned.split() if x], dtype=float)
    except Exception:
        return np.array([], dtype=float)


def extract_protein_id(full_name):
    parts = str(full_name).split('|')
    if len(parts) >= 2:
        return parts[1].strip()
    return str(full_name).split()[0].strip()


def load_disorder_metrics(csv_path, threshold=0.5):
    if not csv_path:
        print(f"  [disorder] SKIP: no path provided")
        return {}
    path = Path(csv_path)
    if not path.exists() or not path.is_file():
        print(f"  [disorder] SKIP: file not found ({csv_path})")
        return {}
    mdf = pd.read_csv(path)
    metrics = {}
    for _, row in mdf.iterrows():
        acc    = extract_protein_id(row["name"])
        scores = parse_disorder_scores(row["disorder_scores"])
        if len(scores) == 0:
            continue
        n_dis_domains = row.get("num_disordered_domains", np.nan)
        n_fol_domains = row.get("num_folded_domains",     np.nan)
        metrics[acc] = {
            "avg_disorder":           float(np.mean(scores)),
            "frac_disordered":        float(np.mean(scores >= threshold)),
            "frac_disordered_length": float(np.sum(scores >= threshold) / len(scores)),
            "num_disordered_domains": float(n_dis_domains) if not pd.isna(n_dis_domains) else np.nan,
            "num_folded_domains":     float(n_fol_domains) if not pd.isna(n_fol_domains) else np.nan,
        }
    return metrics


def merge_disorder_into_df(df, incendi_metrics, ortholog_metrics):
    disorder_params = [
        "avg_disorder", "frac_disordered", "frac_disordered_length",
        "num_disordered_domains", "num_folded_domains",
    ]
    inc_ids  = df["incendi_protein"].apply(extract_protein_id)
    orth_ids = df["ortholog_accession"].apply(extract_protein_id)
    for p in disorder_params:
        df[f"incendi_{p}"]  = inc_ids.map( {k: v[p] for k, v in incendi_metrics.items()})
        df[f"ortholog_{p}"] = orth_ids.map({k: v[p] for k, v in ortholog_metrics.items()})
    return df

def compute_pi_for_ids(fasta_path, target_ids):

    path = Path(fasta_path)
    if not path.exists() or not path.is_file():
        print(f"  [pI] SKIP: file not found ({fasta_path})")
        return {}
    pi_dict = {}
    n_found = 0
    for record in SeqIO.parse(str(path), "fasta"):
        acc = extract_protein_id(record.id)
        if acc not in target_ids:
            continue
        try:
            pi_dict[acc] = IsoelectricPoint(str(record.seq)).pi()
            n_found += 1
        except Exception:
            pass
    return pi_dict


def merge_pi_into_df(df, incendi_pi, ortholog_pi):
    inc_ids  = df["incendi_protein"].apply(extract_protein_id)
    orth_ids = df["ortholog_accession"].apply(extract_protein_id)
    df["incendi_pI"]  = inc_ids.map(incendi_pi)
    df["ortholog_pI"] = orth_ids.map(ortholog_pi)
    return df
# ─── END pI METRICS ───────────────────────────────────────────────────────────

# ─── PARAMS ───────────────────────────────────────────────────────────────────
params = [
    "hbond_per_Rg", "salt_bridges_per_Rg",
    "pipi_per_Rg", "disulfide_per_Rg",
    "ionic_vs_hbond",
    "helix_per_Rg", "strand_per_Rg", "coil_per_Rg",
    "frac_sasa_hydrophobic",
    "frac_sasa_ALA", "frac_sasa_VAL", "frac_sasa_ILE",
    "frac_sasa_LEU", "frac_sasa_MET", "frac_sasa_PRO",
    "frac_sasa_aromatic",
    "frac_sasa_PHE", "frac_sasa_TRP", "frac_sasa_TYR",
    "frac_sasa_nonpolar_uncharged",
    "frac_sasa_polar",
    "frac_sasa_SER", "frac_sasa_THR", "frac_sasa_CYS",
    "frac_sasa_ASN", "frac_sasa_GLN",
    "frac_sasa_charged", 
    "frac_sasa_pos",
    "frac_sasa_LYS", "frac_sasa_ARG",
    "frac_sasa_neg",
    "frac_sasa_ASP", "frac_sasa_GLU",
    "ARG_LYS_ratio", "GLU_ASP_ratio",
    "net_charge",
    "pI",
    "Positive_Surface_Charge_Fraction", "Negative_Surface_Charge_Fraction",
    "frac_disordered", 

    
]

param_labels = {
    "hbond_per_Rg":                     "H-bonds per $R_g$",
    "salt_bridges_per_Rg":              "Salt Bridges per $R_g$",
    "pipi_per_Rg":                      "π–π per $R_g$",
    "disulfide_per_Rg":                 "Disulfides per $R_g$",
    "ionic_vs_hbond":                   "Salt Bridges /H-bonds",
    "helix_per_Rg":                     "Helices per $R_g$",
    "strand_per_Rg":                    "Strands per $R_g$",
    "coil_per_Rg":                      "Coils per $R_g$",
    "frac_sasa_hydrophobic":            "Hydrophobic\nSASA Frac.",
    "frac_sasa_ALA":                    "Ala SASA Frac.",
    "frac_sasa_VAL":                    "Val SASA Frac.",
    "frac_sasa_ILE":                    "Ile SASA Frac.",
    "frac_sasa_LEU":                    "Leu SASA Frac.",
    "frac_sasa_MET":                    "Met SASA Frac.",
    "frac_sasa_PRO":                    "Pro SASA Frac.",
    "frac_sasa_aromatic":               "Aromatic SASA Frac.",
    "frac_sasa_PHE":                    "Phe SASA Frac.",
    "frac_sasa_TRP":                    "Trp SASA Frac.",
    "frac_sasa_TYR":                    "Tyr SASA Frac.",
    "frac_sasa_nonpolar_uncharged":     "Nonpolar SASA Frac.",
    "frac_sasa_polar":                  "Polar SASA Frac.",
    "frac_sasa_SER":                    "Ser SASA Frac.",
    "frac_sasa_THR":                    "Thr SASA Frac.",
    "frac_sasa_CYS":                    "Cys SASA Frac.",
    "frac_sasa_ASN":                    "Asn SASA Frac.",
    "frac_sasa_GLN":                    "Gln SASA Frac.",
    "frac_sasa_charged":                "Charged SASA Frac.",
    "frac_sasa_pos":                    "Pos. SASA Frac.",
    "frac_sasa_LYS":                    "Lys SASA Frac.",
    "frac_sasa_ARG":                    "Arg SASA Frac.",
    "frac_sasa_neg":                    "Neg. SASA Frac.",
    "frac_sasa_ASP":                    "Asp SASA Frac.",
    "frac_sasa_GLU":                    "Glu SASA Frac.",
    "ARG_LYS_ratio":                    "Arg/(Arg+Lys)",
    "GLU_ASP_ratio":                    "Glu/(Glu+Asp)",
    "net_charge":                       "Net Charge",
    "pI":                               "Isoelectric Point",
    "Positive_Surface_Charge_Fraction": "+ Surface Charge Frac.",
    "Negative_Surface_Charge_Fraction": "- Surface Charge Frac.",
    "frac_disordered":                  "Frac. Disordered",
}


incendi_disorder = load_disorder_metrics(INCENDI_DISORDER_CSV, threshold=DISORDER_THRESHOLD)

all_incendi_ids = set()
for csv_path, label, temp, classification in SPECIES:
    path = Path(csv_path)
    if not path.exists():
        continue
    _df = pd.read_csv(csv_path, usecols=["incendi_protein"])
    all_incendi_ids.update(_df["incendi_protein"].apply(extract_protein_id).tolist())

incendi_pi = compute_pi_for_ids(INCENDI_FASTA, all_incendi_ids)

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

    # Disoder scores
    orth_disorder_path = DISORDER_FILES.get(label, "")
    ortholog_disorder  = load_disorder_metrics(orth_disorder_path, threshold=DISORDER_THRESHOLD)
    df = merge_disorder_into_df(df, incendi_disorder, ortholog_disorder)

    #pI
    ortholog_ids    = set(df["ortholog_accession"].apply(extract_protein_id).tolist())
    orth_fasta_path = FASTA_FILES.get(label, "")
    ortholog_pi = compute_pi_for_ids(orth_fasta_path, ortholog_ids)
    df = merge_pi_into_df(df, incendi_pi, ortholog_pi)

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
