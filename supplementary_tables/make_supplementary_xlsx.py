#!/usr/bin/env python3
"""Supplementary Tables S1-S14 of the manuscript as Excel workbooks (Table_S1.xlsx ... Table_S14.xlsx).

Each workbook opens on an 'About' sheet (the table's title, its files and the step of the analysis that made
them), followed by one sheet per file. The values are those of the tab-separated files, unchanged: numbers are
stored as numbers and shown with the same decimals as in the file; everything else is text.

Usage, from the repository's root:
    python3 supplementary_tables/make_supplementary_xlsx.py            # writes supplementary_tables/*.xlsx
    python3 supplementary_tables/make_supplementary_xlsx.py --out DIR
Needs Python 3 with openpyxl.
"""
import argparse
import csv
import re
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

PAPER = ("Long-read haplotypes of the PER3 coding VNTR: chronotype associations track sub-lineages of one "
         "four-unit allele, not repeat length")
REPO = "github.com/comparativechrono/per3-vntr-pangenome"
T = "outputs/tables/"

# number, title (as in the manuscript's list of supplementary material), the step(s) that made the files, files
TABLES = [
    (1, "Data sources and versions, and the script that made each data extract",
     "the scripts in data-build/", ["data-build/data_sources.tsv", "data-build/README.md"]),
    (2, "The repeat in each of the 462 HPRC haplotypes: units by two methods, structure, allele, protein; "
        "agreement with the PGGB graph",
     "steps 02-03", [T + "hprc_haplotype_alleles.tsv", T + "graph_concordance.tsv"]),
    (3, "The 11 unit sequences, their translation and frequency; the eight allele structures",
     "step 02", [T + "unit_catalogue.tsv", T + "allele_structures.tsv"]),
    (4, "Frequencies by population and group, tests of heterogeneity and Hardy-Weinberg equilibrium, "
        "4B sub-lineages",
     "steps 04, 06 and 07", [T + f for f in ["hprc_freq_by_population.tsv", "hprc_lineages_by_group.tsv",
                                             "hprc_tests.tsv", "kg_freq_by_population.tsv",
                                             "kg_lineages_by_population.tsv", "sublineages_by_group.tsv"]]),
    (5, "LD between repeat length and all 5,395 SNVs, by group; the best tags of each allele; extent of LD; "
        "LD between the two tags in 1000 Genomes",
     "steps 05-06", [T + f for f in ["hprc_ld_vntr_snvs.tsv", "hprc_tags_top.tsv", "hprc_lineage_tags.tsv",
                                     "hprc_ld_extent.tsv", "kg_tag_ld.tsv"]]),
    (6, "Validation of tag-based genotypes against long reads in 222 people; discordant people; structural-variant "
        "calls at the repeat in 1000 Genomes, by group and haplotype by haplotype, with the alleles they describe",
     "steps 01, 02 and 06", [T + f for f in ["kg_validation.tsv", "kg_validation_discordant.tsv", "kg_sv_calls.tsv",
                                             "kg_sv_by_group.tsv", "kg_sv_carriers.tsv", "kg_sv_alleles.tsv"]]),
    (7, "The repeat in the phenome-wide scan of Mukamel et al. (2021): sleep traits and the strongest associations "
        "among 787 phenotypes",
     "step 13", [T + "mukamel_sleep_traits.tsv", T + "mukamel_top_phenotypes.tsv"]),
    (8, "Associations of the tags and of the allele markers with 13 traits; accelerometer measures in minutes",
     "step 08", [T + f for f in ["gwas_tag_results.tsv", "gwas_lineage_marker_results.tsv",
                                 "gwas_tag_accelerometer_raw_units.tsv"]]),
    (9, "The alleles and 4B sub-lineages on which the variants of the PER3 chronotype signal lie, in long-read "
        "haplotypes and 1000 Genomes",
     "step 07", [T + f for f in ["lineage_snps_hprc.tsv", "lineage_snps_1kg.tsv", "sublineages_by_group.tsv"]]),
    (10, "Joint models for five traits and three LD panels, with the derived contrasts of the five-unit alleles and "
         "their joint test; stepwise selection by R and GCTA; the tag conditional on selected signals; single-SNP "
         "conditioning; check of the R implementation against GCTA",
     "step 09", [T + f for f in ["lineage_model.tsv", "cojo_selected.tsv", "cojo_tag_conditional.tsv",
                                 "single_snp_conditional.tsv", "cojo_reference_panels.tsv", "gcta_crosscheck.tsv"]]),
    (11, "SuSiE-RSS credible sets, posterior inclusion probabilities of the tags and LD-mismatch diagnostics, with "
         "three LD panels",
     "step 10", [T + f for f in ["susie_credible_sets.tsv", "susie_tags.tsv", "susie_diagnostics.tsv"]]),
    (12, "PER3 eQTL results of GTEx v8 for the tags and the markers of 4B, 4B.1 and 4B.2 in 49 tissues; "
         "associations with plasma PER3 protein in the GWAS Catalog",
     "step 11", [T + f for f in ["gtex_summary.tsv", "gtex_all.tsv", "gwas_catalog_per3_protein.tsv"]]),
    (13, "Units in 35 primate genome sequences with the assembly behind each; units in 53 PER3 proteins; the region "
         "aligned to the human repeat",
     "step 12 (assemblies.tsv: data-build/12_primate_assemblies.py)",
     [T + f for f in ["primate_units_dna.tsv", "primate_units_protein.tsv", "primate_units_by_species.tsv",
                      "primate_aligned_region.tsv"]] + ["data/primates/assemblies.tsv"]),
    (14, "The per-allele difference between five- and four-unit haplotypes implied by the effects of the alleles in "
         "each 1000 Genomes population",
     "step 14", [T + "implied_length_effects.tsv"]),
]
# sheet names longer than Excel's 31 characters, or not unique
SHEET = {"gwas_tag_accelerometer_raw_units.tsv": "gwas_tag_accel_raw_units", "README.md": "data-build_README"}

FONT = Font(name="Arial", size=10)
BOLD = Font(name="Arial", size=10, bold=True)
TITLE = Font(name="Arial", size=14, bold=True)
SUBTITLE = Font(name="Arial", size=11)
HEAD_FILL = PatternFill("solid", fgColor="DDE3EA")
HEAD_LINE = Border(bottom=Side(style="thin", color="7F8C9A"))
NUM = re.compile(r"^-?(0|[1-9]\d*)(?:\.(\d+))?(?:[eE]([-+]?\d+))?$")


def typed(text):
    """A cell's value and number format: numbers as numbers (with the file's decimals), the rest as text."""
    m = NUM.match(text)
    if not m:
        return text, None
    whole, dec, exp = m.groups()
    if dec is None and exp is None:
        return (int(text), "0") if len(whole) <= 15 else (text, None)
    fmt = "0." + "0" * len(dec) if dec else "0"
    return float(text), fmt + ("E+00" if exp is not None else "")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def header_cell(cell):
    cell.font, cell.fill, cell.border = BOLD, HEAD_FILL, HEAD_LINE
    cell.alignment = Alignment(vertical="top")


def add_table_sheet(wb, name, rows):
    ws = wb.create_sheet(name)
    widths = [len(h) for h in rows[0]]
    for j, h in enumerate(rows[0], 1):
        header_cell(ws.cell(row=1, column=j, value=h))
    for i, row in enumerate(rows[1:], 2):
        for j, text in enumerate(row, 1):
            if text == "":
                continue
            value, fmt = typed(text)
            cell = ws.cell(row=i, column=j, value=value)
            if isinstance(value, str):
                cell.data_type = "s"          # never read as a formula
            cell.font = FONT
            if fmt:
                cell.number_format = fmt
            widths[j - 1] = max(widths[j - 1], len(text))
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[ws.cell(row=1, column=j).column_letter].width = min(max(w + 2, 8), 60)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    ws.page_setup.orientation = "landscape"
    ws.print_title_rows = "1:1"             # repeat the header on every printed page
    return len(rows) - 1, len(rows[0])


def add_text_sheet(wb, name, path, lines):
    ws = wb.create_sheet(name)
    header_cell(ws.cell(row=1, column=1, value=path))
    for i, line in enumerate(lines, 2):
        if line:
            cell = ws.cell(row=i, column=1, value=line)
            cell.data_type, cell.font = "s", FONT
    ws.column_dimensions["A"].width = 140
    ws.freeze_panes = "A2"
    return len(lines), 1


def workbook(root, number, title, steps, files):
    wb = Workbook()
    about = wb.active
    about.title = "About"
    contents = []
    for f in files:
        path = root / f
        name = SHEET.get(path.name, path.stem)
        if path.suffix == ".tsv":
            n_rows, n_cols = add_table_sheet(wb, name, read_tsv(path))
        else:
            n_rows, n_cols = add_text_sheet(wb, name, f, path.read_text(encoding="utf-8").splitlines())
        contents.append((name, f, n_rows, n_cols))

    for col, width in zip("ABCD", (24, 110, 9, 9)):
        about.column_dimensions[col].width = width
    about["A1"] = f"Table S{number}."
    about["A1"].font = TITLE
    about.merge_cells("A2:D2")
    about["A2"] = title + "."
    about["A2"].font = SUBTITLE
    about["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    about.row_dimensions[2].height = 15 * (len(title) // 135 + 1) + 3   # rows with wrapped text get a set height
    notes = [("Manuscript", PAPER),
             ("Source", f"{REPO}: tab-separated files made by {steps} of the analysis; values copied unchanged"),
             ("Numbers", "Stored as numbers and shown with the same decimals as in the source files; "
                         "empty cells are missing values")]
    for r, (label, text) in enumerate(notes, 4):
        about.cell(row=r, column=1, value=label).font = BOLD
        cell = about.cell(row=r, column=2, value=text)
        cell.font = FONT
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        about.cell(row=r, column=1).alignment = Alignment(vertical="top")
        about.row_dimensions[r].height = 13 * (len(text) // 115 + 1) + 2
    about.page_setup.orientation = "landscape"
    about.page_setup.fitToWidth, about.page_setup.fitToHeight = 1, 0
    about.sheet_properties.pageSetUpPr.fitToPage = True
    top = 4 + len(notes) + 1
    for j, h in enumerate(["Sheet", "File in the repository", "Rows", "Columns"], 1):
        header_cell(about.cell(row=top, column=j, value=h))
    for i, (name, f, n_rows, n_cols) in enumerate(contents, top + 1):
        for j, v in enumerate([name, f, n_rows, n_cols], 1):
            about.cell(row=i, column=j, value=v).font = FONT
    if any(f.endswith("README.md") for f in files):
        r = top + len(contents) + 2
        cell = about.cell(row=r, column=1, value="The sheet data-build_README holds the text of "
                                                 "data-build/README.md, one line per row.")
        cell.font = FONT
    wb.properties.title = f"Table S{number}"
    wb.properties.creator = "per3-vntr-pangenome"
    return wb


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=here.parent, help="the repository's root")
    ap.add_argument("--out", type=Path, default=here, help="folder for the workbooks")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for number, title, steps, files in TABLES:
        out = args.out / f"Table_S{number}.xlsx"
        workbook(args.root, number, title, steps, files).save(out)
        print(f"{out.name}: {', '.join(Path(f).name for f in files)}")


if __name__ == "__main__":
    main()
