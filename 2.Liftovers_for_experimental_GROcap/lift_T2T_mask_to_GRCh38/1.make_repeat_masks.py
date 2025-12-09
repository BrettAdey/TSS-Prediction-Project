#!/usr/bin/env python3
import sys, re, glob
from Bio import SeqIO

"""
python3 make_repeat_masks.py 'T2T_forward/*.fasta' repeats.bed nonrepeats.bed

Takes soft-masked T2T chromosome fastas as input, one per chromosome.
Outputs .bed files of lowercase (repeat) and uppercase (nonrepeat) runs.
"""

pattern = sys.argv[1]
rep_out = sys.argv[2]
nonrep_out = sys.argv[3]

# Regex patterns for lower and uppercase runs
rep_pat     = re.compile(r'[a-z]+')
nonrep_pat  = re.compile(r'[A-Z]+')

with open(rep_out, 'w') as rep_fh, open(nonrep_out, 'w') as nonrep_fh:
    for fasta in sorted(glob.glob(pattern)):
        for rec in SeqIO.parse(fasta, "fasta"):
            s = str(rec.seq)
            # repeats
            for m in rep_pat.finditer(s):
                rep_fh.write(f"{rec.id}\t{m.start()}\t{m.end()}\n")
            # non-repeats
            for m in nonrep_pat.finditer(s):
                nonrep_fh.write(f"{rec.id}\t{m.start()}\t{m.end()}\n")
print(f"Wrote: {rep_out} and {nonrep_out}")
