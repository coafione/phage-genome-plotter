#!/usr/bin/env python3

import os
import subprocess
import argparse
from itertools import combinations
from pathlib import Path

VALID_EXTENSIONS = [".fasta", ".fa", ".fna"]

def get_phage_ids_from_fasta(fasta_dir):
    phage_ids = []
    for f in os.listdir(fasta_dir):
        if any(f.endswith(ext) for ext in VALID_EXTENSIONS):
            phage_ids.append(Path(f).stem)
    return sorted(set(phage_ids))

def run_blast_all(fasta_dir, output_dir="blast_results", db_dir="blastdb"):
    phage_list = get_phage_ids_from_fasta(fasta_dir)

    if not phage_list:
        print("❌ No valid FASTA files found in:", fasta_dir)
        return

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    print("🔍 Phages detected:", ", ".join(phage_list))
    print("🔬 Creating BLAST databases...")
    for phage in phage_list:
        fasta_path = next(Path(fasta_dir).glob(f"{phage}.*"), None)
        db_path = os.path.join(db_dir, phage)

        if not os.path.exists(f"{db_path}.nin"):
            print(f"→ makeblastdb for {phage}")
            subprocess.run([
                "makeblastdb", "-in", str(fasta_path),
                "-dbtype", "nucl", "-out", db_path
            ])

    print("🚀 Running pairwise BLASTn comparisons...")
    for g1, g2 in combinations(phage_list, 2):
        query = next(Path(fasta_dir).glob(f"{g1}.*"), None)
        db = os.path.join(db_dir, g2)
        out_file = os.path.join(output_dir, f"{g1}_vs_{g2}.txt")

        if not os.path.exists(out_file):
            print(f"BLAST: {g1} → {g2}")
            subprocess.run([
                "blastn", "-query", str(query), "-db", db,
                "-outfmt", "6", "-evalue", "1e-5", "-out", out_file
            ])

    print("✅ All pairwise BLAST jobs completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-run all-vs-all BLASTn for detected phage FASTA files.")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing .fasta/.fa/.fna files")
    parser.add_argument("--output_dir", default="blast_results", help="Directory to store BLAST output")
    parser.add_argument("--db_dir", default="blastdb", help="Directory to store BLAST databases")

    args = parser.parse_args()
    run_blast_all(args.fasta_dir, args.output_dir, args.db_dir)

