import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.stats import false_discovery_control
from pathlib import Path

plt.style.use('default')
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = 14

# ─── CONFIGURATION ────────────────────────────────────────────────────────────
SPECIES = [
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Dictyostelium_discoideum_comparison.csv",  "D. discoideum",    27,  "Meso. Amoeba"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/vvermi_comparison.csv",                    "V. vermiformis",   45,  "Meso./Thermo Tol. Amoeba"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Porphyra_umbilicalis_comparison.csv",      "P. umbilicalis",   34,  "Meso. Algea"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Cyanidioschyzon_merolae_comparison.csv",   "C. merolae",       56,  "Thermo. Algae"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Talaromyces_amestolkiae_comparison.csv",   "T. amestolkiae",   37,  "Meso. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Aspergillus_sergii_comparison.csv",        "A. sergii",        40,  "Meso. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Chaetomium_thermophilum_comparison.csv",   "C. thermophilum",  60,  "Thermo. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Thermomyces_lanuginosus_comparison.csv",   "T. lanuginosus",   63,  "Thermo. Fungi"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Herpetosiphon_gulosus_comparison.csv",     "H. gulosus",       37,  "Meso. Bacteria (Chloroflexi)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Caldilinea_aerophila_comparison.csv",      "C. aerophila",     65,  "Thermo. Bacteria (Chloroflexi)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Deinococcus_ficus_comparison.csv",         "D. ficus",         35,  "Meso. Bacteria (Deinococcota)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Thermus_aquaticus_comparison.csv",         "T. aquaticus",     80,  "Thermo. Bacteria (Deinococcota)"),
    ("/home/jkniblo/Ameoba/struc_params/structural_comparison/Pyrobaculum_ferrireducens_comparison.csv", "P. ferrireducens", 100, "Thermo. Archaea"),
]

OUT_NAME     = "comparison_boxplot_single_parameter_SerSASA_rotated_v2"
TOP_HIT_ONLY = False
# ─── END CONFIGURATION ────────────────────────────────────────────────────────

# ─── DERIVED METRICS ─────────────────────────────────────────────────────────
def add_derived_metrics(df):
    for prefix in ["incendi_", "ortholog_"]:
        p = prefix
        total_surface = df[f"{p}num_exposed_pos"] + df[f"{p}num_exposed_neg"]
        df[f"{p}Positive_Surface_Charge_Fraction"] = df[f"{p}num_exposed_pos"] / total_surface.replace(0, np.nan)
        aro_cols = [f"{p}frac_sasa_TRP", f"{p}frac_sasa_TYR", f"{p}frac_sasa_PHE"]
        present  = [c for c in aro_cols if c in df.columns]
        if present:
            df[f"{p}frac_sasa_aromatic"] = df[present].sum(axis=1)
    return df

# ─── Y-AXIS LABELS ────────────────────────────────────────────────────────────
ylabels = {
    "Positive_Surface_Charge_Fraction": "+ Surface Charge Fraction",
    "frac_sasa_aromatic": "Aromatic SASA Fraction",
    "frac_sasa_SER":      "Serine SASA Fraction",
}

# ─── PARAMS TO PLOT ───────────────────────────────────────────────────────────
params = [
    "frac_sasa_SER",
]

# ─── COLORMAP ────────────────────────────────────────────────────────────────
delta_cmap = LinearSegmentedColormap.from_list(
    'custom',
    ['#7b5ea7', '#fafaf8', '#d4a843'],
    N=256
)
delta_norm = plt.Normalize(vmin=-1, vmax=1)

# ─── LOAD AND PROCESS ALL SPECIES ────────────────────────────────────────────
species_pairs = {param: {} for param in params}
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
    species_labels.append(label)
    species_ns[label]    = len(df)
    species_temp[label]  = temp
    species_class[label] = classification
    print(f"  Loaded {label} ({temp}°C): {len(df)} pairs")

    for param in params:
        col_i = f"incendi_{param}"
        col_o = f"ortholog_{param}"
        if col_i not in df.columns or col_o not in df.columns:
            continue
        iv   = df[col_i].values.astype(float)
        ov   = df[col_o].values.astype(float)
        mask = ~(np.isnan(iv) | np.isnan(ov))
        species_pairs[param][label] = (iv[mask], ov[mask])

params = [
    p for p in params
    if any(len(v[0]) > 0 for v in species_pairs[p].values())
]
print(f"\nPlotting {len(params)} parameters across {len(species_labels)} species")

# ─── FDR: paired t-test (ttest_rel) per (param, species) ─────────────────────
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

# ─── PLOT ─────────────────────────────────────────────────────────────────────
param = params[0]

medians = {}
for label in species_labels:
    pair = species_pairs[param].get(label, (np.array([]), np.array([])))
    iv_p, ov_p = pair
    if len(iv_p) > 0:
        delta = iv_p - ov_p
        medians[label] = np.median(delta)

abs_max = max(abs(v) for v in medians.values()) if medians else 1.0
norm_medians = {label: v / abs_max for label, v in medians.items()}

species_color = {
    label: delta_cmap(delta_norm(norm_medians.get(label, 0.0)))
    for label in species_labels
}

n_sp  = len(species_labels)
box_w = 0.6
rng   = np.random.default_rng(seed=42)

# ── Rotated: wider than tall, species on y-axis ───────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(7, n_sp * 0.37 + 1))
fig.subplots_adjust(left=0.3)

# Species run bottom→top; reverse so first species is at top
y_positions  = np.arange(n_sp)[::-1]
y_ticklabels = []
for label in species_labels:
    italic_label = label.replace(' ', '\\ ')
    y_ticklabels.append(f"$\\it{{{italic_label}}}$")

all_deltas_for_param = []
for label in species_labels:
    pair = species_pairs[param].get(label, (np.array([]), np.array([])))
    iv_p, ov_p = pair
    if len(iv_p) > 0:
        delta = iv_p - ov_p
        all_deltas_for_param.extend(delta)

x_min = min(all_deltas_for_param) if all_deltas_for_param else -1
x_max = max(all_deltas_for_param) if all_deltas_for_param else 1
x_range = x_max - x_min

for s_idx, label in enumerate(species_labels):
    pair = species_pairs[param].get(label, (np.array([]), np.array([])))
    iv_p, ov_p = pair
    if len(iv_p) == 0:
        continue
    
    delta = iv_p - ov_p
    y     = y_positions[s_idx]
    p_fdr = fdr_map[(param, label)]
    col   = species_color[label]

    ax.boxplot(
        delta,
        positions=[y],
        widths=box_w,
        vert=False,                          # ← horizontal boxes
        patch_artist=True,
        medianprops=dict(color='black', linewidth=2),
        boxprops=dict(facecolor=col, alpha=0.85),
        whiskerprops=dict(linewidth=1.5, color='black'),
        capprops=dict(linewidth=1.5, color='black'),
        flierprops=dict(marker=''),
        showfliers=False,
        zorder=2,
    )

    jitter = rng.uniform(-box_w * 0.35, box_w * 0.35, size=len(delta))
    ax.scatter(
        delta, y + jitter,
        color=col, alpha=0.35, s=15,
        edgecolors='k', zorder=1,
    )

    # ── significance star to the right of the whisker ──────────────────────
    if p_fdr < 0.01:
        delta_max = np.max(delta)
        star_x = delta_max + 0.08 * x_range
        ax.text(star_x, y, '★', ha='center', va='center',
                fontsize=18, color='black', zorder=100)

ax.tick_params(which='minor', length=0)
ax.tick_params(which='major', length=6, width=1.5)
ax.axvline(0, color='black', linewidth=1.2, linestyle='--', alpha=0.5)
ax.set_yticks(y_positions)
ax.set_yticklabels(y_ticklabels, fontsize=12)
ax.set_xlabel(f"Δ {ylabels.get(param, param)}", fontsize=12)
ax.set_ylim(-0.6, n_sp - 0.4)
ax.set_xlim(-0.3, x_max + 0.15 * x_range)
ax.spines[["top", "right"]].set_visible(False)

# ─── COLORBAR ─────────────────────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap=delta_cmap, norm=delta_norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.08)
cbar.set_ticks([-1, 0, 1])
cbar.set_ticklabels(['Higher in\nOrtholog', 'No\ndifference', 'Higher in\nIncendi'], fontsize=12)
cbar.ax.tick_params(labelsize=14)

plt.tight_layout()
plt.savefig(f"{OUT_NAME}.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUT_NAME}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUT_NAME}.svg", dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()
print(f"Saved: {OUT_NAME}.png / .svg")