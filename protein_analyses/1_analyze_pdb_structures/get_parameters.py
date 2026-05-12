import sys
import os
import csv
import shutil
import traceback
import pandas as pd
from pathlib import Path
import struc_params as sp
import mdtraj as md
import MDAnalysis as mda

Species = "Thermomyces_lanuginosus"

FILTERED_PAIRS   = Path(f"/home/jkniblo/Ameoba/foldseek/{Species}/filtered_pairs.tsv")
INCENDI_STRUCTS  = Path("/home/jkniblo/Ameoba/AF_predicted/incendi_strucs")
ORTHOLOG_STRUCTS = Path(f"/home/jkniblo/Ameoba/AF_predicted/{Species}_strucs")
TMP_DIR          = Path("/home/jkniblo/Ameoba/tmp_hydrogens")
OUT_CSV          = f"{Species}_comparison.csv"
TMSCORE_MIN      = 0.5


TMP_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)


# ─── HELPERS ─────────────────────────────────────────────────────────────────

def analyze_pdb(pdb_path: Path, tmp_dir: Path) -> dict:
    
    original_dir = os.getcwd()
    os.chdir(tmp_dir)

    try:
        pdb_file = str(pdb_path)
        mdA_u = mda.Universe(pdb_file)
        mdT_u = md.load_pdb(pdb_file)

        num_chains, num_residues, rg, helix_num, strand_num, coil_num = \
            sp.analyze_geometry(mdT_u)

        total_sasa, hydro_frac, polar_frac, charged_frac, frac_pos, frac_neg, frac_sasa_per_residue = \
            sp.analyze_surface(mdT_u)

        num_hbonds, num_salt_bridge, num_pipi, num_disulfides = \
            sp.analyze_interactions(pdb_file, mdA_u)

        net_charge, surface_positive, surface_negative, surface_net_charge, highly_exposed_frac = \
            sp.analyze_electrostats(mdT_u)

    except Exception as e:
        print(f"    ERROR analyzing {pdb_path.name}: {e}")
        traceback.print_exc()
        return None
    finally:
        os.chdir(original_dir)

    result = {
        "num_chains":           num_chains,
        "num_residues":         num_residues,
        "rg_nm":                float(rg),
        "helix_num":            int(helix_num),
        "strand_num":           int(strand_num),
        "coil_num":             int(coil_num),
        "total_sasa_A2":        float(total_sasa),
        "frac_sasa_hydrophobic": float(hydro_frac),
        "frac_sasa_polar":      float(polar_frac),
        "frac_sasa_charged":    float(charged_frac),
        "frac_sasa_pos":        float(frac_pos),
        "frac_sasa_neg":        float(frac_neg),
        "num_hbonds":           int(num_hbonds),
        "num_salt_bridges":     int(num_salt_bridge),
        "num_disulfides":       int(num_disulfides),
        "num_pipi":             int(num_pipi),
        "net_charge":           int(net_charge),
        "num_exposed_pos":      int(surface_positive),
        "num_exposed_neg":      int(surface_negative),
        "surface_net_charge":   int(surface_net_charge),
        "highly_exposed_frac":  float(highly_exposed_frac),
    }

    for res, frac in frac_sasa_per_residue.items():
        result[f"frac_sasa_{res}"] = float(frac[0]) if hasattr(frac, '__len__') else float(frac)

    return result


def prefix_dict(d: dict, prefix: str) -> dict:
    return {f"{prefix}{k}": v for k, v in d.items()}



pairs = []
with open(FILTERED_PAIRS) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        tmscore = float(row["alntmscore"])
        if tmscore < TMSCORE_MIN:
            continue
        pairs.append({
            "incendi_protein":    row["incendi_protein"].strip(),
            "ortholog_accession": row["ortholog_accession"].strip(),
            "alntmscore":         tmscore,
            "pident":             float(row["pident"]),
        })


rows = []
skipped = 0

for i, pair in enumerate(pairs):
    incendi_id = pair["incendi_protein"]
    target_id = pair["ortholog_accession"]

    incendi_pdb = INCENDI_STRUCTS  / f"{incendi_id}.pdb"
    target_pdb = ORTHOLOG_STRUCTS / f"{target_id}.pdb"

    if not incendi_pdb.exists():
        print(f"  SKIP: incendi PDB not found: {incendi_pdb}")
        skipped += 1
        continue
    if not target_pdb.exists():
        print(f"  SKIP: ortholog PDB not found: {target_pdb}")
        skipped += 1
        continue

    print(f"  Analyzing incendi:  {incendi_pdb.name}")
    incendi_results = analyze_pdb(incendi_pdb, TMP_DIR)
    if incendi_results is None:
        skipped += 1
        continue

    print(f"  Analyzing ortholog: {target_pdb.name}")
    target_results = analyze_pdb(target_pdb, TMP_DIR)
    if target_results is None:
        skipped += 1
        continue

    row = {
        "incendi_protein":    incendi_id,
        "ortholog_accession": target_id,
        "alntmscore":         pair["alntmscore"],
        "pident":             pair["pident"],
        **prefix_dict(incendi_results, "incendi_"),
        **prefix_dict(target_results, "ortholog_"),
    }
    rows.append(row)


df = pd.DataFrame(rows)

# add difference columns for key features (incendi - ortholog)
diff_cols = [
    "num_residues", "rg_nm", "helix_num", "strand_num", "coil_num",
    "total_sasa_A2", "frac_sasa_hydrophobic", "frac_sasa_polar",
    "frac_sasa_charged", "num_hbonds", "num_salt_bridges",
    "num_disulfides", "num_pipi", "net_charge", "surface_net_charge",
    "highly_exposed_frac",
]
for col in diff_cols:
    if f"incendi_{col}" in df.columns and f"ortholog_{col}" in df.columns:
        df[f"diff_{col}"] = df[f"incendi_{col}"] - df[f"ortholog_{col}"]

df.to_csv(OUT_CSV, index=False)
shutil.rmtree(TMP_DIR, ignore_errors=True)
