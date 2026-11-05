# Supplementary tables as Excel workbooks

`Table_S1.xlsx` to `Table_S14.xlsx` are Supplementary Tables S1–S14 of the manuscript, one workbook per table.
Each opens on an **About** sheet – the table's title, its files and the step of the analysis that made them –
followed by one sheet per file, in the order the manuscript lists them.

The values are those of the tab-separated files in `outputs/tables/`, `data-build/` and `data/primates/`,
unchanged: numbers are stored as numbers and shown with the same decimals as in the files, and everything else is
text. Empty cells are missing values. `sublineages_by_group.tsv` belongs to Tables S4 and S9 and is in both, and
Table S1 has the text of `data-build/README.md` on a sheet of its own, one line per row. The tab-separated files
remain the primary record.

To remake the workbooks after `run_all.sh` (needs Python 3 with openpyxl):

```bash
python3 supplementary_tables/make_supplementary_xlsx.py
```
