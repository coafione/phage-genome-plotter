#!/usr/bin/env python3

import json
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.lines import Line2D

def identity_to_color(identity, identity_min=80, colormap="viridis"):
    cmap = cm.get_cmap(colormap)
    norm = Normalize(vmin=identity_min, vmax=100)
    return cmap(norm(identity))

def plot_genomes(data, phage_order, output_file="phage_comparison.svg", output_format="svg",
                 cds_height=1.6, gap=4.5, min_cds_length=50, blast_identity_min=80,
                 colormap="viridis", fig_width=18, fig_height=14, bar_thickness=6, font_size=10):

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_xlim(0, max(d["length"] for d in data.values()) + 1000)
    ax.set_ylim(-gap * len(phage_order), 2)
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax.set_xlabel("Genomic position (bp)", fontsize=font_size)

    phage_pos = {phage: -i * gap for i, phage in enumerate(phage_order)}

    # Draw genomes and CDS
    for phage in phage_order:
        y = phage_pos[phage]
        length = data[phage]["length"]

        ax.plot([0, length], [y, y], lw=bar_thickness, color="gray")
        ax.text(length + 700, y, phage, va='center', ha='left', fontsize=font_size)

        for start, end, strand, _ in data[phage]["cds"]:
            if end - start < min_cds_length:
                continue
            color = "green" if strand == 1 else "orange"
            box_y = y + cds_height * 0.6
            rect = patches.Rectangle(
                (start, box_y - cds_height / 2), end - start, cds_height,
                linewidth=0.5, edgecolor='black', facecolor=color, alpha=0.85
            )
            ax.add_patch(rect)

    # Draw BLAST links
    for phage in phage_order:
        src_y = phage_pos[phage]
        for link in data[phage]["links"]:
            target = link["target"]
            if target not in phage_pos:
                continue
            tgt_y = phage_pos[target]
            src_idx = phage_order.index(phage)
            tgt_idx = phage_order.index(target)

            if abs(src_idx - tgt_idx) != 1:
                continue
            if src_idx > tgt_idx:
                continue

            poly = patches.Polygon([
                [link["qstart"], src_y - 0.4],
                [link["qend"], src_y - 0.4],
                [link["send"], tgt_y + 0.4],
                [link["sstart"], tgt_y + 0.4],
            ], closed=True, facecolor=identity_to_color(link["identity"], blast_identity_min, colormap), alpha=0.5)
            ax.add_patch(poly)

    # Colorbar
    cax = fig.add_axes([0.93, 0.35, 0.015, 0.4])
    norm = Normalize(vmin=blast_identity_min, vmax=100)
    cb = ColorbarBase(cax, cmap=cm.get_cmap(colormap), norm=norm, orientation='vertical')
    cb.set_label(f'BLAST Identity (%)\n≥ {blast_identity_min}', fontsize=font_size)

    # Strand color legend
    legend_elements = [
        Line2D([0], [0], color='green', lw=6, label='+ Strand CDS'),
        Line2D([0], [0], color='orange', lw=6, label='– Strand CDS')
    ]
    fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.015),
               ncol=2, fontsize=font_size-1, frameon=True)

    fig.suptitle("Comparative Genomics", fontsize=font_size+6, y=0.97)
    plt.tight_layout(rect=[0, 0.05, 0.9, 0.95])
    plt.savefig(output_file, format=output_format, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot comparative genomics of phages based on parsed data.")
    parser.add_argument("--input_json", required=True, help="Input JSON file from parser")
    parser.add_argument("--output_file", default="phage_comparison.svg", help="Output file name (default: SVG)")
    parser.add_argument("--format", default="svg", help="Output file format: svg, png, pdf, etc.")
    parser.add_argument("--phage_order", default=None, help="Comma-separated list of phages for plot order (optional)")

    # New flexible arguments
    parser.add_argument("--colormap", default="viridis", help="Colormap for BLAST links (default: viridis)")
    parser.add_argument("--fig_width", type=float, default=18, help="Figure width in inches (default: 18)")
    parser.add_argument("--fig_height", type=float, default=14, help="Figure height in inches (default: 14)")
    parser.add_argument("--bar_thickness", type=float, default=6, help="Thickness of genome bars (default: 6)")
    parser.add_argument("--font_size", type=int, default=10, help="Font size for labels (default: 10)")

    args = parser.parse_args()

    with open(args.input_json) as f:
        data = json.load(f)

    if args.phage_order:
        phage_order = [p.strip() for p in args.phage_order.split(",")]
    else:
        phage_order = list(data.keys())

    plot_genomes(data, phage_order,
                 args.output_file, args.format,
                 cds_height=1.6, gap=4.5, min_cds_length=50, blast_identity_min=80,
                 colormap=args.colormap, fig_width=args.fig_width, fig_height=args.fig_height,
                 bar_thickness=args.bar_thickness, font_size=args.font_size)
