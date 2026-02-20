#!/usr/bin/env python

import subprocess
import os
import pandas as pd

metadata = pd.read_csv('sample_metadata.csv')
metadata['name'] = metadata['brain_region'] + "_" + metadata['stress_condition'] + "_" + metadata['replicate']
rename_dict = dict(zip(metadata["run"], metadata["name"]))

SCRATCH_SRA = "/tscc/lustre/ddn/scratch/b2long/bulkRNAseq_APA/sra"

SCRATCH_FASTQ = "/tscc/lustre/ddn/scratch/b2long/bulkRNAseq_APA/raw_fastq"
os.makedirs(SCRATCH_FASTQ, exist_ok=True)

for sra_id, sample_name in rename_dict.items():
    print("Currently downloading: " + sra_id)
    subprocess.call(f"prefetch {sra_id}", shell=True)

    print("Generating fastq for: " + sra_id)
    subprocess.call(f"fasterq-dump --outdir {SCRATCH_FASTQ} {sra_id}", shell=True)

    # rename and gzip any file containing the SRA ID
    for fname in os.listdir(SCRATCH_FASTQ):
        if sra_id in fname and fname.endswith(".fastq"):
            new_fname = fname.replace(sra_id, sample_name)
            src = os.path.join(SCRATCH_FASTQ, fname)
            dst = os.path.join(SCRATCH_FASTQ, new_fname)
            os.rename(src, dst)
            print(f"Renamed {fname} -> {new_fname}")
            subprocess.call(f"gzip {dst}", shell=True)
            print(f"Gzipped {new_fname}")
    
    # delete the prefetched directory to save space
    sra_file = os.path.join(SCRATCH_SRA, f"{sra_id}.sra")
    if os.path.exists(sra_file):
        os.remove(sra_file)
        print(f"Deleted prefetch cache for: {sra_id}")
    else:
        print(f"No prefetch directory found for: {sra_id}, skipping cleanup")