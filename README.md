# Phage Genome Plotter 🧬

**Phage Genome Plotter** is a flexible, customizable toolkit for comparative genomics visualization.  
It combines **nucleotide-level similarity** (via BLASTn) with **annotation data** (from GenBank .gbk files), allowing users to see not just *where* genomes align, but also *what* biological features (genes, proteins) are involved.
---

## ✨ Features
- Runs all-vs-all `blastn` from FASTA files
- Parses GenBank `.gbk` + BLAST results automatically
- Produces beautiful, scalable genome plots
- Flexible control over colors, figure size, and labels
- Simple command-line interface (CLI)
- Outputs to `.svg`, `.png`, `.pdf`, etc.

---

## 🚀 Quickstart

1. **Install Requirements**
```bash
pip install biopython matplotlib
sudo apt install ncbi-blast+  # (or conda install -c bioconda blast)




## 📚 Citation

If you use **Phage Genome Plotter** for your research, please cite it as:

> Constanza Afione Di Cristofano, *Phage Genome Plotter: a customizable toolkit for comparative genomics visualization*, GitHub (2025). [https://github.com/coafione/phage-genome-plotter](https://github.com/yourusername/phage-genome-plotter)

Thank you for supporting open science! 🧬✨

