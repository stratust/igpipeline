#!/usr/bin/env python
import os
import sys
import glob
from Bio import SeqIO
import re

data_dir, out_fname = sys.argv[1:]
data_dir = os.path.normpath(data_dir)
out_fname = os.path.normpath(out_fname)

ab1_files = glob.glob(data_dir+"/*.ab1")
print("- converting", len(ab1_files), "ab1 files to fastq ...")

seqs = []
for f in ab1_files:
    f = os.path.normpath(f)
    record = SeqIO.read(f, "abi")
    # use filename instead of abi internal names (otherwise names will be
    # different if files are renamed)
    record.id = re.search('(.*)\.ab1', os.path.basename(f)).groups()[0]
    seqs.append(record)

count = SeqIO.write(seqs, out_fname, "fastq-illumina")

assert len(ab1_files) == count, "conversion not successful for all sequences"

