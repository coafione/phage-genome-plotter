#!/usr/bin/env python3

import os
import argparse
import json
from pathlib import Path
from Bio import SeqIO

def parse_gbk_files(gbk_dir, valid_ids):
    data = {}
    valid_gbk_exts = [".gbk", ".gbff", ".gb"]
    for file in os.listdir(gbk_dir):
        if not any(file.endswith(ext) for ext in valid_gbk_exts):
            continue
        phage_id = Path(file).stem
        if phage_id not in valid_ids:
            continue

        path = os.path.join(gbk_dir, file)
        record = SeqIO.read(path, "genbank")
        cds_list = []

        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                gene = feature.qualifiers.get("gene", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                locus = feature.qualifiers.get("locus_tag", [""])[0]
                label = gene or product or locus
                cds_list.append((start, end, strand, label))

        data[phage_id] = {
            "length": len(record.seq),
            "cds": cds_list,
            "links": []
        }

    return data

def parse_blast_files(blast_dir, data, identity_thresh=80.0, min_len=100):
    for file in os.listdir(blast_dir):
        if not file.endswith(".txt"):
            continue
        parts = file.replace(".txt", "").split("_vs_")
        if len(parts) != 2:
            continue
        g1, g2 = parts
        if g1 not in data or g2 not in data:
            continue

        with open(os.path.join(blast_dir, file)) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 10:
                    continue

                identity = float(cols[2])
                align_len = abs(int(cols[7]) - int(cols[6]))
                if identity < identity_thresh or align_len < min_len:
                    continue

                qstart, qend = sorted([int(cols[6]), int(cols[7])])
                sstart, send = sorted([int(cols[8]), int(cols[9])])
                data[g1]["links"].append({
                    "target": g2,
                    "qstart": qstart,
                    "qend": qend,
                    "sstart": sstart,
                    "send": send,
                    "identity": identity
                })

def auto_detect_phage_ids(fasta_dir):
    ids = []
    for f in os.listdir(fasta_dir):
        if any(f.endswith(ext) for ext in [".fasta", ".fa", ".fna"]):
            ids.append(Path(f).stem)
    return sorted(set(ids))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse .gbk and BLAST results into a unified JSON for plotting.")
    parser.add_argument("--gbk_dir", required=True, help="Directory containing GenBank .gbk files")
    parser.add_argument("--blast_dir", required=True, help="Directory with BLAST .txt files (outfmt 6)")
    parser.add_argument("--fasta_dir", required=True, help="Directory with .fasta/.fa/.fna files (used for ID detection)")
    parser.add_argument("--output_file", default="parsed_data.json", help="Output JSON file")

    args = parser.parse_args()
    phage_ids = auto_detect_phage_ids(args.fasta_dir)
    data = parse_gbk_files(args.gbk_dir, phage_ids)
    parse_blast_files(args.blast_dir, data)
    
    with open(args.output_file, "w") as f:
        json.dump(data, f, indent=2)

    print(f"✅ Parsed data saved to {args.output_file}")

