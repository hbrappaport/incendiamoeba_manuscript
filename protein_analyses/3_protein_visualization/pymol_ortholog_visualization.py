"""
Publication-quality PyMOL grid: incendi protein + orthologs, single row.
Designed for Cell journal figure standards (300 DPI, ray-traced, PDF output).
"""

import io
import os
import numpy as np
import pandas as pd
from pathlib import Path
from pymol import cmd, stored
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas as rl_canvas
from reportlab.lib.utils import ImageReader

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

TARGET_PROTEIN = "g19137.t2"   # <- SET THIS

AF_BASE     = "/home/jkniblo/Ameoba/AF_predicted"
INCENDI_DIR = f"{AF_BASE}/incendi_strucs"
CSV_DIR     = "/home/jkniblo/Ameoba/struc_params/structural_comparison"
OUT_DIR     = "nicer_images_v2"
TEMP_DIR    = "/tmp/pymol_pub_frames"

# (tsv_filename, display_name, pdb_subfolder)
SPECIES = [
    ("Dictyostelium_discoideum_comparison.csv",  "D. discoideum",    "Dictyostelium_discoideum_strucs"),
    ("vvermi_comparison.csv",                    "V. vermiformis",   "vvermi_strucs"),
    ("Porphyra_umbilicalis_comparison.csv",      "P. umbilicalis",   "Porphyra_umbilicalis_strucs"),
    ("Cyanidioschyzon_merolae_comparison.csv",   "C. merolae",       "Cyanidioschyzon_merolae_strucs"),
    ("Talaromyces_amestolkiae_comparison.csv",   "T. amestolkiae",   "Talaromyces_amestolkiae_strucs"),
    ("Aspergillus_sergii_comparison.csv",        "A. sergii",        "Aspergillus_sergii_strucs"),
    ("Chaetomium_thermophilum_comparison.csv",   "C. thermophilum",  "Chaetomium_thermophilum_strucs"),
    ("Herpetosiphon_gulosus_comparison.csv",     "H. gulosus",       "Herpetosiphon_gulosus_strucs"),
    ("Caldilinea_aerophila_comparison.csv",      "C. aerophila",     "Caldilinea_aerophila_strucs"),
    ("Deinococcus_ficus_comparison.csv",         "D. ficus",         "Deinococcus_ficus_strucs"),
    ("Thermus_aquaticus_comparison.csv",         "T. aquaticus",     "Thermus_aquaticus_strucs"),
    ("Pyrobaculum_ferrireducens_comparison.csv", "P. ferrireducens", "Pyrobaculum_ferrireducens_strucs"),
]

# ── Color helpers ─────────────────────────────────────────────────────────────
def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def hex_to_pymol(h):
    return "0x" + h.lstrip("#").upper()

# ── Charge colors — edit the three CSS hex strings only ──────────────────────
POS_HEX_CSS = "#3A6AA1"   # basic residues  (Lys, Arg)
NEG_HEX_CSS = "#C1532D"   # acidic residues (Asp, Glu)
NEU_HEX_CSS = "#FFFFFF"   # neutral

POS_RGB = hex_to_rgb(POS_HEX_CSS)
NEG_RGB = hex_to_rgb(NEG_HEX_CSS)
NEU_RGB = hex_to_rgb(NEU_HEX_CSS)

POS_HEX = hex_to_pymol(POS_HEX_CSS)
NEG_HEX = hex_to_pymol(NEG_HEX_CSS)
NEU_HEX = hex_to_pymol(NEU_HEX_CSS)

# ── pLDDT discrete bands — edit hex strings to remap colors ──────────────────
# (b_low, b_high_exclusive, hex_css)
PLDDT_BANDS = [
    (0,   50,  "#B97878"),   # Very low  — dusty rose
    (50,  70,  "#C8AF78"),   # Low       — warm sand
    (70,  90,  "#6E9BC3"),   # Confident — steel blue
    (90,  101, "#416AA0"),   # Very high — slate blue
]

# ── Ortholog selection ────────────────────────────────────────────────────────
TM_THRESHOLD = 0.5   # minimum alntmscore; set to 0.0 to disable

# ── Layout ────────────────────────────────────────────────────────────────────
OUTPUT_DPI = 300
CELL_PX    = 900
LABEL_H    = 80
PANEL_GAP  = 0
MARGIN     = 30

# ── Typography ────────────────────────────────────────────────────────────────
FONT_PATHS = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
]
FONT_BOLD_PATHS = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
]
FONT_SIZE_SPECIES = 22
FONT_SIZE_ACC     = 18
FONT_SIZE_META    = 17
FONT_SIZE_TITLE   = 24

# ── Canvas colors ─────────────────────────────────────────────────────────────
BG             = (255, 255, 255)
TEXT_PRIMARY   = (25,  25,  25)
TEXT_SECONDARY = (90,  90,  90)

RENDER_MODES = [
    ("cartoon", "plddt",  "cartoon_plddt",  ""),
    ("cartoon", "charge", "cartoon_charge", ""),
    # ("surface", "plddt",  "surface_plddt",  ""),
    # ("surface", "charge", "surface_charge", ""),
]

# Shared camera view — locked from incendi, reused for all panels
SHARED_VIEW = None

# ══════════════════════════════════════════════════════════════════════════════
# UTILITY
# ══════════════════════════════════════════════════════════════════════════════

def load_font(paths, size):
    for p in paths:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


# ══════════════════════════════════════════════════════════════════════════════
# ORTHOLOG GATHERING
# Reads pident (0–100) and alntmscore directly from comparison TSV.
# ══════════════════════════════════════════════════════════════════════════════

def gather_orthologs(pid: str) -> list:
    """
    Per species: keep hits with alntmscore > TM_THRESHOLD, then select
    the single ortholog with the lowest pident (most distant hit).
    """
    orthologs = []
    for csv_fname, sp_name, sp_dir in SPECIES:
        csv_path = Path(CSV_DIR) / csv_fname
        if not csv_path.exists():
            print(f"  SKIP {sp_name}: file not found ({csv_path})")
            continue

        try:
            df = pd.read_csv(str(csv_path), sep=None, engine="python")
        except Exception as e:
            print(f"  SKIP {sp_name}: could not read {csv_path}: {e}")
            continue

        if "ortholog_accession" not in df.columns:
            print(f"  SKIP {sp_name}: missing 'ortholog_accession' column")
            continue

        df  = df.dropna(subset=["ortholog_accession"])
        sub = df[df["incendi_protein"].astype(str).str.strip() == pid]
        if sub.empty:
            continue

        candidates = []
        for _, row in sub.iterrows():
            acc = str(row["ortholog_accession"]).strip()
            pdb = Path(AF_BASE) / sp_dir / f"{acc}.pdb"

            tm = None
            if "alntmscore" in row:
                try:
                    tm = float(row["alntmscore"])
                except Exception:
                    pass

            try:
                fid = float(row["pident"])   # already 0–100
            except Exception:
                print(f"    Could not read pident for {acc}, skipping")
                continue

            if tm is not None and tm < TM_THRESHOLD:
                continue

            candidates.append({
                "sp_name":  sp_name,
                "acc":      acc,
                "pdb_path": str(pdb),
                "homology": fid,
                "tmscore":  tm,
            })

        if not candidates:
            print(f"  No hits passing TM>{TM_THRESHOLD} for {sp_name}")
            continue

        best   = min(candidates, key=lambda x: x["homology"])
        tm_str = f"{best['tmscore']:.3f}" if best["tmscore"] is not None else "n/a"
        print(f"  {sp_name:22s}  {best['acc']:20s}  "
              f"pident={best['homology']:.1f}%  TM={tm_str}")
        orthologs.append(best)

    return orthologs


# ══════════════════════════════════════════════════════════════════════════════
# SHARED VIEW
# ══════════════════════════════════════════════════════════════════════════════

def compute_shared_view(incendi_pdb: str) -> tuple:
    """Zoom to incendi, capture camera, clean up."""
    cmd.reinitialize()
    cmd.load(incendi_pdb, "__ref__")
    cmd.zoom("__ref__", buffer=4)
    view = cmd.get_view()
    cmd.reinitialize()
    return view


# ══════════════════════════════════════════════════════════════════════════════
# COLORING — uses alter() to set color indices at atom level.
# This bypasses PyMOL's ribbon interpolation entirely.
# ══════════════════════════════════════════════════════════════════════════════

def _register_color(name: str, hex_css: str) -> int:
    """Register a named color and return its PyMOL index."""
    cmd.set_color(name, [c / 255.0 for c in hex_to_rgb(hex_css)])
    return cmd.get_color_index(name)


def apply_charge(obj: str):
    neu_idx = _register_color("_neu_c", NEU_HEX_CSS)
    pos_idx = _register_color("_pos_c", POS_HEX_CSS)
    neg_idx = _register_color("_neg_c", NEG_HEX_CSS)

    # Set neutral on everything first, then override charged residues
    cmd.alter(obj, f"color = {neu_idx}")
    cmd.alter(f"({obj}) and resn LYS+ARG", f"color = {pos_idx}")
    cmd.alter(f"({obj}) and resn ASP+GLU", f"color = {neg_idx}")
    cmd.rebuild()


def apply_plddt(obj: str):
    # Set the lowest-confidence color as default (catches b==0 edge case)
    default_idx = _register_color("_plddt_default", PLDDT_BANDS[0][2])
    cmd.alter(obj, f"color = {default_idx}")

    for i, (lo, hi, css) in enumerate(PLDDT_BANDS):
        idx = _register_color(f"_plddt_{i}", css)
        # PyMOL selection algebra uses > not >=
        cmd.alter(f"({obj}) and b > {lo - 1} and b < {hi}", f"color = {idx}")

    cmd.rebuild()


# ══════════════════════════════════════════════════════════════════════════════
# PYMOL RENDERING
# ══════════════════════════════════════════════════════════════════════════════

def _pymol_publication_base():
    cmd.set("ray_trace_mode",        1)
    cmd.set("ray_shadows",           1)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias",             2)
    cmd.set("ray_trace_gain",        0.1)
    cmd.set("ambient",               0.45)
    cmd.set("direct",                0.65)
    cmd.set("reflect",               0.20)
    cmd.set("shininess",             25)
    cmd.set("specular",              0.20)
    cmd.set("depth_cue",             1)
    cmd.set("fog_start",             0.55)
    cmd.bg_color("white")


def _apply_cartoon_settings(discrete=False):
    cmd.set("cartoon_discrete_colors",     1)   
    cmd.set("cartoon_smooth_loops",        0)
    cmd.set("cartoon_fancy_helices",       0)   # back on
    cmd.set("cartoon_cylindrical_helices", 0)   # back off
    cmd.set("cartoon_loop_radius",         0.25)
    cmd.set("cartoon_tube_radius",         0.25)
    cmd.set("cartoon_oval_length",         1.4)
    cmd.set("cartoon_rect_length",         1.4)
    cmd.set("cartoon_rect_width",          0.35)
    cmd.set("cartoon_helix_radius",        1.1)


def _apply_surface_settings():
    cmd.set("ray_trace_mode",  1)    # Phong — realistic, not cel-shaded
    cmd.set("surface_quality", 3)
    cmd.set("surface_solvent", 0)
    cmd.set("transparency",    0.0)
    cmd.set("surface_mode",    3)    # per-atom coloring
    cmd.set("ray_trace_gain",  0.05)
    cmd.set("specular",        0.10)
    cmd.set("shininess",       10)
    cmd.set("ambient",         0.55)
    cmd.set("direct",          0.60)




def render_structure(pdb_path: str, safe_name: str,
                     representation: str, coloring: str,
                     ref_pdb: str = None):
    tmp = Path(TEMP_DIR) / f"{safe_name}_{representation}_{coloring}.png"
    if tmp.exists():
        return str(tmp)

    cmd.reinitialize()
    try:
        cmd.load(pdb_path, safe_name)
    except Exception as e:
        print(f"    Cannot load {pdb_path}: {e}")
        return None

    # Align ortholog to incendi reference
    if ref_pdb:
        ref = "__incendi_ref__"
        try:
            cmd.load(ref_pdb, ref)
            try:
                cmd.align(safe_name, ref, quiet=1)
            except Exception:
                try:
                    cmd.super(safe_name, ref, quiet=1)
                except Exception:
                    pass
            cmd.delete(ref)
        except Exception:
            pass

    _pymol_publication_base()
    cmd.hide("everything", "all")

    discrete = coloring in ("charge", "plddt")

    if representation == "cartoon":
        # Apply settings BEFORE show so geometry is built correctly
        _apply_cartoon_settings(discrete=discrete)
        cmd.show("cartoon", safe_name)
    elif representation == "surface":
        cmd.show("surface", safe_name)
        _apply_surface_settings()

    # Color via alter() — atom-level, bypasses ribbon interpolation
    if coloring == "plddt":
        apply_plddt(safe_name)
    elif coloring == "charge":
        apply_charge(safe_name)

    # Lock to shared view so all panels are framed identically
    if SHARED_VIEW is not None:
        cmd.set_view(SHARED_VIEW)
    else:
        cmd.zoom(safe_name, buffer=4)

    cmd.png(str(tmp), width=CELL_PX, height=CELL_PX, dpi=OUTPUT_DPI, ray=1)
    cmd.reinitialize()

    return str(tmp) if tmp.exists() else None


# ══════════════════════════════════════════════════════════════════════════════
# GRID COMPOSITOR -> PDF
# ══════════════════════════════════════════════════════════════════════════════

def build_and_save_pdf(panels: list, pdf_path: str,
                       coloring: str, mode_title: str):
    fnt_sp  = load_font(FONT_BOLD_PATHS, FONT_SIZE_SPECIES)
    fnt_acc = load_font(FONT_PATHS,      FONT_SIZE_ACC)
    fnt_met = load_font(FONT_PATHS,      FONT_SIZE_META)
    fnt_ttl = load_font(FONT_BOLD_PATHS, FONT_SIZE_TITLE)

    n       = len(panels)
    panel_h = CELL_PX + LABEL_H

    total_w = MARGIN + n * CELL_PX + MARGIN
    total_h = MARGIN + panel_h + MARGIN

    fig  = Image.new("RGB", (total_w, total_h), BG)
    draw = ImageDraw.Draw(fig)

    x_cursor = MARGIN

    for idx, p in enumerate(panels):
        y = MARGIN

        # Structure image
        if p["img_path"] and Path(p["img_path"]).exists():
            try:
                struc_img = (Image.open(p["img_path"])
                                  .convert("RGB")
                                  .resize((CELL_PX, CELL_PX), Image.LANCZOS))
                fig.paste(struc_img, (x_cursor, y))
            except Exception as e:
                print(f"    Paste error {p['img_path']}: {e}")
                draw.rectangle([x_cursor, y,
                                x_cursor + CELL_PX, y + CELL_PX],
                               fill=(220, 220, 220))
        else:
            draw.rectangle([x_cursor, y, x_cursor + CELL_PX, y + CELL_PX],
                           fill=(220, 220, 220))
            draw.text((x_cursor + 10, y + CELL_PX // 2),
                      "render failed", font=fnt_acc, fill=(140, 140, 140))

        # Labels on white background — no grey bar
        strip_y = y + CELL_PX

        # Species name
        bb = draw.textbbox((0, 0), p["sp_name"], font=fnt_sp)
        draw.text(
            (x_cursor + (CELL_PX - (bb[2] - bb[0])) // 2, strip_y + 8),
            p["sp_name"], font=fnt_sp, fill=TEXT_PRIMARY)

        # pident — same size as species name
        if p.get("homology") is not None:
            met = f"{p['homology']:.1f}% identity"
            bb  = draw.textbbox((0, 0), met, font=fnt_sp)
            draw.text(
                (x_cursor + (CELL_PX - (bb[2] - bb[0])) // 2,
                 strip_y + 8 + FONT_SIZE_SPECIES + 4),
                met, font=fnt_sp, fill=TEXT_PRIMARY)

        x_cursor += CELL_PX

    # Mode title
    draw.text((MARGIN, 6), mode_title, font=fnt_ttl, fill=TEXT_SECONDARY)

    # Write PDF via reportlab (lossless TIFF embedded)
    px_to_pt  = 72.0 / OUTPUT_DPI
    page_w_pt = total_w * px_to_pt
    page_h_pt = total_h * px_to_pt

    pdf_canvas = rl_canvas.Canvas(pdf_path, pagesize=(page_w_pt, page_h_pt))

    buf = io.BytesIO()
    fig.save(buf, format="TIFF", compression="tiff_lzw",
             dpi=(OUTPUT_DPI, OUTPUT_DPI))
    buf.seek(0)

    pdf_canvas.drawImage(ImageReader(buf), 0, 0,
                         width=page_w_pt, height=page_h_pt,
                         preserveAspectRatio=False)

    pdf_canvas.setTitle(Path(pdf_path).stem)
    pdf_canvas.setAuthor("Sukenik Lab")
    pdf_canvas.setSubject("AlphaFold ortholog structural comparison")
    pdf_canvas.setCreator("visualize_orthologs_pub.py / PyMOL + reportlab")
    pdf_canvas.save()

    print(f"  -> {Path(pdf_path).name}  "
          f"({n} panels | {page_w_pt:.0f} x {page_h_pt:.0f} pt | {OUTPUT_DPI} DPI)")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUT_DIR,  exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)

pid = TARGET_PROTEIN
print(f"\nTarget protein: {pid}")

incendi_pdb = Path(INCENDI_DIR) / f"{pid}.pdb"
if not incendi_pdb.exists():
    raise FileNotFoundError(f"Incendi PDB not found: {incendi_pdb}")

# Lock shared camera view from incendi structure before any rendering
SHARED_VIEW = compute_shared_view(str(incendi_pdb))
print(f"Shared view locked from {pid}\n")

print("Gathering orthologs:")
orthologs = gather_orthologs(pid)
if not orthologs:
    raise RuntimeError(f"No orthologs found for {pid}. Check CSVs and species list.")

print(f"\nSelected {len(orthologs)} ortholog(s) after TM>{TM_THRESHOLD} filter")

structures = [
    {"sp_name":    "I. cascadensis",
     "acc":        pid,
     "pdb_path":   str(incendi_pdb),
     "homology":   None,
     "tmscore":    None,
     "ref_pdb":    None,
     "is_incendi": True}
]
for o in orthologs:
    if Path(o["pdb_path"]).exists():
        structures.append({**o, "ref_pdb": str(incendi_pdb), "is_incendi": False})
    else:
        print(f"  Missing PDB skipped: {o['pdb_path']}")

print(f"\nRendering {len(structures)} structures x 4 modes = "
      f"{len(structures) * 4} ray-traced panels\n")

for representation, coloring, suffix, mode_title in RENDER_MODES:
    print(f"Mode: {mode_title}")
    panels = []

    for s in structures:
        safe = s["acc"][:18].replace(" ", "_")
        img  = render_structure(
            pdb_path       = s["pdb_path"],
            safe_name      = safe,
            representation = representation,
            coloring       = coloring,
            ref_pdb        = s.get("ref_pdb"),
        )
        status = "OK" if img else "FAIL"
        print(f"  [{status}] {s['sp_name']:22s} {s['acc'][:22]}")
        panels.append({**s, "img_path": img})

    out_pdf = str(Path(OUT_DIR) / f"{pid}_{suffix}.pdf")
    build_and_save_pdf(panels, out_pdf, coloring, mode_title)

print(f"\nDone. 4 PDFs saved -> {OUT_DIR}")