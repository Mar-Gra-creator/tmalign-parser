#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Script for serial TMalign (single-threaded version)

Version: 1.2 
Author: Marcin Gradowski
Contact: marcin_gradowski@sggw.edu.pl

Requirements:
  - TM-align in PATH
  - Python 3
  - pandas, regex, glob, os, re, csv, datetime

Usage:
  1. Place query PDBs in 'q/' and template PDBs in 't/'
  2. Run: python tmalign_parser.py

Outputs:
  - c-tm/, c-tm/all, c-tm/all_tab, c-tm/all_tab_sorted, tm-stacks/
  - summary tables with filenames only (e.g., e1bamA1.pdb, e2vldA2.pdb)
'''
import os
import glob
import time
import pandas as pd
import re
from datetime import datetime

# Analysis parameters
analysis_name = "my_job"
date_str = datetime.now().strftime('%Y%m%d')
score_threshold = 0.6  # TM-score cutoff for stacks


def create_folder(path):
    os.makedirs(path, exist_ok=True)


def timer(start, end):
    h, rem = divmod(end - start, 3600)
    m, s = divmod(rem, 60)
    return f"{int(h):02}:{int(m):02}:{s:05.2f}"


def run_tm_align(query, template_dir, out_dir):
    """Run TM-align of query against all templates."""
    create_folder(out_dir)
    for templ in glob.glob(os.path.join(template_dir, '*.pdb')):
        qname = os.path.basename(query)
        tname = os.path.basename(templ)
        outfile = os.path.join(out_dir, f"{os.path.splitext(qname)[0]}_{os.path.splitext(tname)[0]}.aln")
        os.system(f"TMalign {query} {templ} > {outfile}")


def parse_and_save():
    """Parse .aln files, extract metrics, sort, filter, and save results with filenames only."""
    # Prepare directories
    create_folder('c-tm/all_tab')
    create_folder('c-tm/all_tab_sorted')
    create_folder('tm-stacks')

    all_summary = []
    stack_counts = {}

    for alnfile in glob.glob('c-tm/all/*.aln'):
        basename = os.path.splitext(os.path.basename(alnfile))[0]
        tab_out = f'c-tm/all_tab/{basename}.tab'
        with open(alnfile) as f:
            text = f.read()
        blocks = re.split(r'\*+\s*TM-align.*?\*+', text)

        records = []
        for block in blocks[1:]:
            # Extract only filenames (no directories)
            # Query filename from basename
            query_file = f"{basename}.pdb"
            # Template filename from parsed ID
            t_m = re.search(r'Name of Chain_2:\s*([^\s(]+)', block)
            template_id = t_m.group(1) if t_m else ''
            template_file = f"{template_id}.pdb"

            # Extract metrics
            qlen_m = re.search(r'Length of Chain_1:\s*(\d+)', block)
            tlen_m = re.search(r'Length of Chain_2:\s*(\d+)', block)
            rmsd_m = re.search(r'RMSD=\s*([\d\.]+)', block)
            id_m = re.search(r'Seq_ID=n_identical/n_aligned=\s*([\d\.]+)', block)
            s1_m = re.search(r'TM-score=\s*([\d\.]+)\s*\(if normalized by length of Chain_1,.*?\)', block)
            s2_m = re.search(r'TM-score=\s*([\d\.]+)\s*\(if normalized by length of Chain_2,.*?\)', block)

            qlen = qlen_m.group(1) if qlen_m else ''
            tlen = tlen_m.group(1) if tlen_m else ''
            RMSD = rmsd_m.group(1) if rmsd_m else ''
            ident = id_m.group(1) if id_m else ''
            score1 = s1_m.group(1) if s1_m else ''
            score2 = s2_m.group(1) if s2_m else ''
            try:
                chosen = max(float(score1) if score1 else 0.0, float(score2) if score2 else 0.0)
            except:
                chosen = 0.0

            records.append([query_file, template_file, qlen, tlen, RMSD, ident, score1, score2, chosen])

        df = pd.DataFrame(records, columns=['query','template','qlen','tlen','RMSD','ident','TM1','TM2','chosen'])
        df.to_csv(tab_out, sep='\t', index=False)

        # Sort and save
        df['chosen'] = pd.to_numeric(df['chosen'], errors='coerce')
        df_sorted = df.sort_values('chosen', ascending=False)
        sorted_out = f'c-tm/all_tab_sorted/{basename}_sorted.tab'
        create_folder(os.path.dirname(sorted_out))
        df_sorted.to_csv(sorted_out, sep='\t', index=False)
        all_summary.append(df_sorted)

        # Filter stacks
        stack = df_sorted[df_sorted['chosen'] >= score_threshold]
        stack_counts[basename] = len(stack)
        if not stack.empty:
            stack.to_csv(f'tm-stacks/{basename}.tab', sep='\t', index=False, header=False)

    # Combined outputs
    combo = pd.concat(all_summary, ignore_index=True)
    combo.to_csv(f'{analysis_name}_all_results_{date_str}.tab', sep='\t', index=False)
    sc = pd.DataFrame(list(stack_counts.items()), columns=['query','count'])
    sc.to_csv(f'{analysis_name}_stack_counts_{date_str}.tab', sep='\t', index=False)


def main():
    # Run TM-align serially
    create_folder('c-tm/all')
    for query in glob.glob('q/*.pdb'):
        run_tm_align(query, 't', f'c-tm/{os.path.splitext(os.path.basename(query))[0]}')
        os.system(f"cat c-tm/{os.path.splitext(os.path.basename(query))[0]}/*.aln > c-tm/all/{os.path.splitext(os.path.basename(query))[0]}.aln")
    parse_and_save()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f"Done in {timer(start_time, time.time())}")

