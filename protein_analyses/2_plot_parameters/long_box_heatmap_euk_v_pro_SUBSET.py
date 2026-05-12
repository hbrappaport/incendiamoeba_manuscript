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

# ─── CONFIGURATION ────────────────────────────────────────────────────────────
SPECIES = [
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Dictyostelium_discoideum_comparison.csv",  "D. discoideum",    27,  "Meso. Amoeba"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/vvermi_comparison.csv",                    "V. vermiformis",   45,  "Meso./Thermo Tol. Amoeba"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Porphyra_umbilicalis_comparison.csv",      "P. umbilicalis",   34,  "Meso. Algea"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Cyanidioschyzon_merolae_comparison.csv",   "C. merolae",       56,  "Thermo. Algae"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Talaromyces_amestolkiae_comparison.csv",   "T. amestolkiae",   37,  "Meso. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Aspergillus_sergii_comparison.csv",        "A. sergii",        40,  "Meso. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Chaetomium_thermophilum_comparison.csv",   "C. thermophilum / T. thermophila",  60,  "Thermo. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Thermomyces_lanuginosus_comparison.csv",   "T. lanuginosus",   63,  "Thermo. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Herpetosiphon_gulosus_comparison.csv",     "H. gulosus",       37,  "Meso. Bacteria (Chloroflexi)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Caldilinea_aerophila_comparison.csv",      "C. aerophila",     65,  "Thermo. Bacteria (Chloroflexi)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Deinococcus_ficus_comparison.csv",         "D. ficus",         35,  "Meso. Bacteria (Deinococcota)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Thermus_aquaticus_comparison.csv",         "T. aquaticus",     80,  "Thermo. Bacteria (Deinococcota)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Pyrobaculum_ferrireducens_comparison.csv", "P. ferrireducens", 100, "Thermo. Archaea"),
]

DOMAIN_SPLIT_AFTER = 7
OUT_NAME           = "delta_heatmap_by_species_by_euk_pro_withIDR_v2"
TOP_HIT_ONLY       = False

# ─── DISORDER CONFIGURATION ───────────────────────────────────────────────────
INCENDI_DISORDER_CSV = "/home/jkniblo/Ameoba/quantify_disorder/incendi_500_metapredictv3_Mar2026.csv"

DISORDER_FILES = {
    "D. discoideum":    "/home/jkniblo/Ameoba/quantify_disorder/Dictyostelium_discoideum_metapredictv3_Mar2026.csv",
    "V. vermiformis":   "/home/jkniblo/Ameoba/quantify_disorder/vvermi_metapredictv3_Mar2026.csv",
    "P. umbilicalis":   "/home/jkniblo/Ameoba/quantify_disorder/Porphyra_umbilicalis_metapredictv3_Mar2026.csv",
    "C. merolae":       "/home/jkniblo/Ameoba/quantify_disorder/Cyanidioschyzon_merolae_metapredictv3_Mar2026.csv",
    "T. amestolkiae":   "/home/jkniblo/Ameoba/quantify_disorder/Talaromyces_amestolkiae_metapredictv3_Mar2026.csv",
    "A. sergii":        "/home/jkniblo/Ameoba/quantify_disorder/Aspergillus_sergii_metapredictv3_Mar2026.csv",
    "C. thermophilum":  "/home/jkniblo/Ameoba/quantify_disorder/Chaetomium_thermophilum_metapredictv3_Mar2026.csv",
    "T. lanuginosus":   "/home/jkniblo/Ameoba/quantify_disorder/Thermomyces_lanuginosus_metapredictv3_Mar2026.csv",
    "H. gulosus":       "/home/jkniblo/Ameoba/quantify_disorder/Herpetosiphon_gulosus_metapredictv3_Mar2026.csv",
    "C. aerophila":     "/home/jkniblo/Ameoba/quantify_disorder/Caldilinea_aerophila_metapredictv3_Mar2026.csv",
    "D. ficus":         "/home/jkniblo/Ameoba/quantify_disorder/Deinococcus_ficus_metapredictv3_Mar2026.csv",
    "T. aquaticus":     "/home/jkniblo/Ameoba/quantify_disorder/Thermus_aquaticus_metapredictv3_Mar2026.csv",
    "P. ferrireducens": "/home/jkniblo/Ameoba/quantify_disorder/Pyrobaculum_ferrireducens_metapredictv3_Mar2026.csv",
}
DISORDER_THRESHOLD = 0.5   # standard metapredict convention
# ─── END CONFIGURATION ────────────────────────────────────────────────────────

nonpolar_uncharged = ["ALA", "VAL", "ILE", "LEU", "MET"]

# ─── DERIVED METRICS ──────────────────────────────────────────────────────────
def add_derived_metrics(df):
    for prefix in ["incendi_", "ortholog_"]:
        p = prefix
        df[f"{p}residue_per_Rg"]      = df[f"{p}num_residues"]     / df[f"{p}rg_nm"]
        df[f"{p}sasa_per_residue"]    = df[f"{p}total_sasa_A2"]    / df[f"{p}num_residues"]
        df[f"{p}hbond_per_Rg"]        = df[f"{p}num_hbonds"]       / df[f"{p}rg_nm"]
        df[f"{p}salt_bridges_per_Rg"] = df[f"{p}num_salt_bridges"] / df[f"{p}rg_nm"]
        df[f"{p}pipi_per_Rg"]         = df[f"{p}num_pipi"]         / df[f"{p}rg_nm"]
        df[f"{p}disulfide_per_Rg"]    = df[f"{p}num_disulfides"]   / df[f"{p}rg_nm"]
        df[f"{p}helix_per_Rg"]        = df[f"{p}helix_num"]        / df[f"{p}rg_nm"]
        df[f"{p}strand_per_Rg"]       = df[f"{p}strand_num"]       / df[f"{p}rg_nm"]
        df[f"{p}coil_per_Rg"]         = df[f"{p}coil_num"]         / df[f"{p}rg_nm"]
        df[f"{p}helix_per_residue"]   = df[f"{p}helix_num"]        / df[f"{p}num_residues"]
        df[f"{p}strand_per_residue"]  = df[f"{p}strand_num"]       / df[f"{p}num_residues"]
        df[f"{p}coil_per_residue"]    = df[f"{p}coil_num"]         / df[f"{p}num_residues"]
        df[f"{p}hbonds_per_helix"]        = df[f"{p}num_hbonds"]       / df[f"{p}helix_num"].replace(0, np.nan)
        df[f"{p}hbonds_per_strand"]       = df[f"{p}num_hbonds"]       / df[f"{p}strand_num"].replace(0, np.nan)
        df[f"{p}salt_bridges_per_helix"]  = df[f"{p}num_salt_bridges"] / df[f"{p}helix_num"].replace(0, np.nan)
        df[f"{p}salt_bridges_per_strand"] = df[f"{p}num_salt_bridges"] / df[f"{p}strand_num"].replace(0, np.nan)
        df[f"{p}pipi_per_helix"]          = df[f"{p}num_pipi"]         / df[f"{p}helix_num"].replace(0, np.nan)
        df[f"{p}pipi_per_strand"]         = df[f"{p}num_pipi"]         / df[f"{p}strand_num"].replace(0, np.nan)
        df[f"{p}disulfide_per_SS"]        = df[f"{p}num_disulfides"] / (
            df[f"{p}helix_num"] + df[f"{p}strand_num"] + df[f"{p}coil_num"]
        )
        npu_cols = [f"{p}frac_sasa_{aa}" for aa in nonpolar_uncharged if f"{p}frac_sasa_{aa}" in df.columns]
        df[f"{p}frac_sasa_nonpolar_uncharged"] = df[npu_cols].sum(axis=1)
        df[f"{p}ionic_vs_hbond"] = df[f"{p}num_salt_bridges"] / df[f"{p}num_hbonds"].replace(0, np.nan)
        total_surface = df[f"{p}num_exposed_pos"] + df[f"{p}num_exposed_neg"]
        df[f"{p}surface_charge_density"]           = total_surface / df[f"{p}num_residues"]
        df[f"{p}surface_charge_per_SASA"]          = total_surface / df[f"{p}total_sasa_A2"]
        df[f"{p}Positive_Surface_Charge_Fraction"] = df[f"{p}num_exposed_pos"] / total_surface.replace(0, np.nan)
        df[f"{p}Negative_Surface_Charge_Fraction"] = df[f"{p}num_exposed_neg"] / total_surface.replace(0, np.nan)
        df[f"{p}Net_Surface_Charge"]               = df[f"{p}surface_net_charge"]
        df[f"{p}ARG_LYS_ratio"]  = df[f"{p}frac_sasa_ARG"] / (df[f"{p}frac_sasa_ARG"] + df[f"{p}frac_sasa_LYS"]).replace(0, np.nan)
        df[f"{p}GLU_ASP_ratio"]  = df[f"{p}frac_sasa_GLU"] / (df[f"{p}frac_sasa_GLU"] + df[f"{p}frac_sasa_ASP"]).replace(0, np.nan)
        df[f"{p}surface_pos_vs_total_charge"] = df[f"{p}num_exposed_pos"] / (abs(df[f"{p}net_charge"]) + 1e-6)
        df[f"{p}surface_neg_vs_total_charge"] = df[f"{p}num_exposed_neg"] / (abs(df[f"{p}net_charge"]) + 1e-6)
        df[f"{p}surface_charge_vs_total"] = (df[f"{p}num_exposed_pos"] + df[f"{p}num_exposed_neg"]) / (abs(df[f"{p}net_charge"]) + 1e-6)
        aro_cols = [f"{p}frac_sasa_TRP", f"{p}frac_sasa_TYR", f"{p}frac_sasa_PHE"]
        present  = [c for c in aro_cols if c in df.columns]
        if present:
            df[f"{p}frac_sasa_aromatic"] = df[present].sum(axis=1)
    return df

# ─── DISORDER METRICS ─────────────────────────────────────────────────────────
def parse_disorder_scores(score_str):
    """Parse metapredict disorder_scores string → 1D numpy array."""
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
    """
    Pull the bare accession from a full FASTA header, e.g.
      'tr|I0HYG1|I0HYG1_CALAS Uncharacterized...'  →  'I0HYG1'
    Falls back to the first whitespace-delimited token if no pipes are found.
    """
    parts = str(full_name).split('|')
    if len(parts) >= 2:
        return parts[1].strip()
    return str(full_name).split()[0].strip()


def load_disorder_metrics(csv_path, threshold=0.5):
    """
    Load a metapredict CSV and return a dict:
        accession → {avg_disorder, frac_disordered,
                     num_disordered_domains, num_folded_domains,
                     frac_disordered_length}
    """
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
        "avg_disorder",
        "frac_disordered",
        "frac_disordered_length",
        "num_disordered_domains",
        "num_folded_domains",
    ]
    inc_ids  = df["incendi_protein"].apply(extract_protein_id)
    orth_ids = df["ortholog_accession"].apply(extract_protein_id)
    for p in disorder_params:
        df[f"incendi_{p}"]  = inc_ids.map( {k: v[p] for k, v in incendi_metrics.items()})
        df[f"ortholog_{p}"] = orth_ids.map({k: v[p] for k, v in ortholog_metrics.items()})
    return df
# ─── END DISORDER METRICS ─────────────────────────────────────────────────────

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
    "rg_nm":                            "Radius of\nGyration",
    "total_sasa_A2":                    "Total\nSASA",
    "helix_num":                        "Number of\nHelices",
    "strand_num":                       "Number of\nStrands",
    "coil_num":                         "Number of\nCoils",
    "num_hbonds":                       "H-bonds",
    "num_salt_bridges":                 "Salt\nBridges",
    "num_disulfides":                   "Disulfides",
    "num_pipi":                         "π–π\nInteractions",
    "net_charge":                       "Net Charge",
    "num_exposed_pos":                  "Exposed\nPos. Residues",
    "num_exposed_neg":                  "Exposed\nNeg. Residues",
    "surface_net_charge":               "Surface\nNet Charge",
    "highly_exposed_frac":              "Highly Exposed\nResidue Frac.",
    "frac_sasa_hydrophobic":            "Hydrophobic SASA Frac.",
    "frac_sasa_polar":                  "Polar SASA Frac.",
    "frac_sasa_charged":                "Charged\nSASA Frac.",
    "frac_sasa_pos":                    "Positive\nSASA Frac.",
    "frac_sasa_neg":                    "Negative\nSASA Frac.",
    "frac_sasa_GLY":                    "Glycine\nSASA Frac.",
    "frac_sasa_ALA":                    "Alanine\nSASA Frac.",
    "frac_sasa_VAL":                    "Valine\nSASA Frac.",
    "frac_sasa_ILE":                    "Isoleucine\nSASA Frac.",
    "frac_sasa_LEU":                    "Leucine\nSASA Frac.",
    "frac_sasa_PRO":                    "Proline\nSASA Frac.",
    "frac_sasa_PHE":                    "Phenylalanine\nSASA Frac.",
    "frac_sasa_TRP":                    "Tryptophan\nSASA Frac.",
    "frac_sasa_TYR":                    "Tyrosine\nSASA Frac.",
    "frac_sasa_MET":                    "Methionine\nSASA Frac.",
    "frac_sasa_SER":                    "Serine\nSASA Frac.",
    "frac_sasa_THR":                    "Threonine\nSASA Frac.",
    "frac_sasa_CYS":                    "Cysteine\nSASA Frac.",
    "frac_sasa_ASN":                    "Asparagine\nSASA Frac.",
    "frac_sasa_GLN":                    "Glutamine\nSASA Frac.",
    "frac_sasa_ASP":                    "Aspartate SASA Frac.",
    "frac_sasa_GLU":                    "Glutamate SASA Frac.",
    "frac_sasa_LYS":                    "Lysine SASA Frac.",
    
    "frac_sasa_ARG":                    "Arginine SASA Frac.",
    "frac_sasa_HIS":                    "Histidine\nSASA Frac.",
    "residue_per_Rg":                   "Residues per\nRadius of Gyration",
    "sasa_per_residue":                 "SASA per\nResidue",
    "hbond_per_Rg":                     "H-bonds per $R_g$",
    "salt_bridges_per_Rg":              "Salt Bridges per $R_g$",
    "pipi_per_Rg":                      "π–π per\nRadius of Gyration",
    "disulfide_per_Rg":                 "Disulfides per\nRadius of Gyration",
    "helix_per_Rg":                     "Helices per\nRadius of Gyration",
    "strand_per_Rg":                    "Strands per\nRadius of Gyration",
    "coil_per_Rg":                      "Coils per\nRadius of Gyration",
    "helix_per_residue":                "Helix\nFraction",
    "strand_per_residue":               "Strand\nFraction",
    "coil_per_residue":                 "Coil\nFraction",
    "hbonds_per_helix":                 "H-bonds\nper Helix",
    "hbonds_per_strand":                "H-bonds\nper Strand",
    "salt_bridges_per_helix":           "Salt Bridges\nper Helix",
    "salt_bridges_per_strand":          "Salt Bridges\nper Strand",
    "pipi_per_helix":                   "π–π\nper Helix",
    "pipi_per_strand":                  "π–π\nper Strand",
    "disulfide_per_SS":                 "Disulfides per\nSecondary Structure",
    "ionic_vs_hbond":                   "Salt Bridges /\nH-bonds Ratio",
    "surface_charge_density":           "Surface Charge\nDensity",
    "surface_charge_per_SASA":          "Surface Charge\nper SASA",
    "Net_Surface_Charge":               "Net Surface\nCharge",
    "Positive_Surface_Charge_Fraction": '+ Surface Charge Frac.', #"Positive Surface\nCharge Frac.",
    "Negative_Surface_Charge_Fraction": '- Surface Charge Frac.',#"Negative Surface\nCharge Frac.",
    "frac_sasa_aromatic":               "Aromatic SASA Frac.",
    "frac_sasa_nonpolar_uncharged":     "Nonpolar SASA Frac.",
    "avg_disorder":                     "Mean\nDisorder Score",
    "frac_disordered":                  "Fraction\nDisordered",
    "frac_disordered_length":           "Disordered\nLength Frac.",
    "num_disordered_domains":           "No. Disordered\nDomains",
    "num_folded_domains":               "No. Folded\nDomains",
    "frac_sasa_SER":                    "Serine SASA Frac.",

}

# ─── LOAD DATA ────────────────────────────────────────────────────────────────
print("Loading incendi disorder metrics...")
incendi_disorder = load_disorder_metrics(INCENDI_DISORDER_CSV, threshold=DISORDER_THRESHOLD)
print(f"  {len(incendi_disorder)} incendi proteins with disorder data")

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
    if TOP_HIT_ONLY:
        df = (df.sort_values("alntmscore", ascending=False)
                .drop_duplicates("incendi_protein")
                .reset_index(drop=True))

    df = add_derived_metrics(df)

    # ── merge disorder scores ──────────────────────────────────────────────────
    orth_disorder_path = DISORDER_FILES.get(label, "")
    ortholog_disorder  = load_disorder_metrics(orth_disorder_path, threshold=DISORDER_THRESHOLD)
    print(f"  [disorder] {label}: {len(ortholog_disorder)} ortholog proteins loaded")
    df = merge_disorder_into_df(df, incendi_disorder, ortholog_disorder)

    species_labels.append(label)
    species_ns[label]    = len(df)
    species_temp[label]  = temp
    species_class[label] = classification
    print(f"  Loaded {label} ({temp}°C): {len(df)} pairs")

    # ── 1:1 ortholog count ────────────────────────────────────────────────────
    counts   = df.groupby("incendi_protein")["ortholog_accession"].count()
    n_single = (counts == 1).sum()
    n_unique = len(counts)
    print(f"    1:1 pairs: {n_single}/{n_unique} ({100*n_single/n_unique:.1f}%) "
          f"of matched incendi proteins")

    # ── homology summary ──────────────────────────────────────────────────────
    pident_vals = df["pident"].dropna()
    print(f"    Avg. % identity: {pident_vals.mean():.2f}%  "
          f"(median: {pident_vals.median():.2f}%,  "
          f"range: {pident_vals.min():.1f}–{pident_vals.max():.1f}%,  "
          f"n={len(pident_vals)})")

    # ── store paired arrays for paired t-test ─────────────────────────────────
    for param in params:
        col_i = f"incendi_{param}"
        col_o = f"ortholog_{param}"
        if col_i not in df.columns or col_o not in df.columns:
            continue
        iv   = df[col_i].values.astype(float)
        ov   = df[col_o].values.astype(float)
        mask = ~(np.isnan(iv) | np.isnan(ov))
        species_pairs[param][label] = (iv[mask], ov[mask])

# ─── OVERALL HOMOLOGY SUMMARY ─────────────────────────────────────────────────
print("\n── One-to-one pairs per species ──")
all_pident   = []
all_single   = 0
all_proteins = 0
for csv_path, label, temp, classification in SPECIES:
    path = Path(csv_path)
    if not path.exists():
        continue
    _df = pd.read_csv(csv_path, usecols=["pident", "incendi_protein", "ortholog_accession"])
    all_pident.append(_df["pident"].dropna())
    counts       = _df.groupby("incendi_protein")["ortholog_accession"].count()
    single_count = (counts == 1).sum()
    total_count  = len(counts)
    all_single   += single_count
    all_proteins += total_count
    print(f"  {label:25s}: {single_count:4d} / {total_count:4d} one-to-one  "
          f"({100*single_count/total_count:5.1f}%)")

all_pident = pd.concat(all_pident)
print(f"\n── Overall statistics across all species ──")
print(f"  Mean % identity:   {all_pident.mean():.2f}%")
print(f"  Median % identity: {all_pident.median():.2f}%")
print(f"  Range:             {all_pident.min():.1f}–{all_pident.max():.1f}%")
print(f"  Total ortholog pairs (rows): {len(all_pident)}")
print(f"  Total unique incendi proteins (across species): {all_proteins}")
print(f"  Proteins with exactly 1 ortholog: {all_single}/{all_proteins} "
      f"({100*all_single/all_proteins:.1f}%)")

# Drop params with no data in any species
params = [
    p for p in params
    if any(len(v[0]) > 0 for v in species_pairs[p].values())
]
print(f"\n{len(params)} parameters, {len(species_labels)} species")

# ─── FDR — paired t-test (ttest_rel) ─────────────────────────────────────────
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

# ─── BUILD MATRICES ───────────────────────────────────────────────────────────
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

# ─── PER-COLUMN NORMALIZATION TO [-1, 1] ──────────────────────────────────────
norm_matrix = np.full_like(median_matrix, np.nan)
for j in range(n_params):
    col = median_matrix[:, j]
    col_abs_max = np.nanmax(np.abs(col))
    if col_abs_max > 0:
        norm_matrix[:, j] = col / col_abs_max
    else:
        norm_matrix[:, j] = 0.0

# ─── TEMPERATURE COLORMAP for y-axis labels ───────────────────────────────────
all_temps = [species_temp[l] for l in species_labels]
temp_norm = plt.Normalize(vmin=min(all_temps), vmax=max(all_temps))
temp_cmap = plt.cm.bwr_r

# ─── PLOT ─────────────────────────────────────────────────────────────────────
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

# ── significance stars ────────────────────────────────────────────────────────
for i in range(n_species):
    for j in range(n_params):
        if sig_matrix[i, j]:
            ax.text(j, i, '★', ha='center', va='center',
                    fontsize=16, color='black', zorder=3)

# ── x-axis ────────────────────────────────────────────────────────────────────
x_labels = [param_labels.get(p, p) for p in params]
ax.set_xticks(range(n_params))
ax.set_xticklabels(x_labels, fontsize=16, rotation=90)

# ── y-axis ────────────────────────────────────────────────────────────────────
ax.set_yticks(range(n_species))
y_labels = []
for label in species_labels:
    italic_label   = label.replace(' ', '\\ ')
    classification = species_class[label]
    y_labels.append(f"$\\it{{{italic_label}}}$\n({classification}, {species_temp[label]}°C)")
ax.set_yticklabels(y_labels, fontsize=16)
for ytick, label in zip(ax.get_yticklabels(), species_labels):
    ytick.set_color('k')

# ── grid lines ────────────────────────────────────────────────────────────────
ax.set_xticks(np.arange(n_params) - 0.5, minor=True)
ax.set_yticks(np.arange(n_species) - 0.5, minor=True)
ax.grid(which='minor', color='white', linewidth=0.8)
ax.tick_params(which='minor', length=0)
ax.tick_params(which='major', length=6, width=1.5)

# ── eukaryote / prokaryote domain divider ─────────────────────────────────────
if DOMAIN_SPLIT_AFTER is not None:
    ax.axhline(
        y=DOMAIN_SPLIT_AFTER + 0.5,
        color='black',
        linewidth=2.0,
        linestyle='--',
        zorder=4,
    )

# ── colorbar ──────────────────────────────────────────────────────────────────
cbar = fig.colorbar(im, ax=ax, pad=0.01)
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels(
    ['Higher in\nOrtholog', 'No\ndifference', 'Higher in\nIncendi'],
    fontsize=14,
)
cbar.ax.tick_params(length=6, width=1.5)

fig.tight_layout()
plt.savefig(f"{OUT_NAME}.png", dpi=200, bbox_inches='tight')
plt.savefig(f"{OUT_NAME}.svg", dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.savefig(f"{OUT_NAME}.pdf", dpi=200, bbox_inches='tight', pad_inches=0.1)

plt.close()
print(f"Saved: {OUT_NAME}.png / .svg / .pdf")