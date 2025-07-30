# tmalign_parser.py

**Version:** 1.2

A lightweight script for serial TM-align of PDB files. Compares each file in `q/` against all in `t/`, parses TM-scores and other metrics, and saves results as TSV.

## Requirements

- Python 3.10
- TM-align in `PATH`
- Python packages: pandas, regex, glob, os, re, csv, datetime

## Usage

1. Place query PDBs in `q/` and template PDBs in `t/`.
2. Run:
   ```bash
   python tmalign_parser.py
   ```

## Output

- `c-tm/all/` – combined `.aln` files per query
- `c-tm/all_tab/` – parsed TSV files
- `c-tm/all_tab_sorted/` – TSV files sorted by TM-score
- `tm-stacks/` – lines with TM-score ≥ threshold (default 0.4)
- Summary tables: `<analysis_name>_all_results_<date>.tab`, `<analysis_name>_stack_counts_<date>.tab`

## Configuration

Edit the variables at the top of the script:
- `analysis_name` (output prefix)
- `score_threshold` (TM-score cutoff)

## License

BSD 2-Clause "Simplified"
