import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
import argparse
import os

def plot_coverage_with_annotations(file_path, pathogen="sars-cov-2"):
    # Load data
    df = pd.read_csv(file_path, sep="\t")

    # Compute bin midpoint and label
    df["bin_mid"] = (df["start"] + df["end"]) // 2
    df = df.sort_values("bin_mid")
    df["bin_label"] = df["bin_mid"].astype(str)

    # Sample type (assuming "EM" or "research" in sample names)
    df["type"] = df["sample"].apply(lambda x: "EM" if "EM" in x else "research")

    # Gene annotations by pathogen
    if pathogen.lower() == "sars-cov-2":
        gene_annotations = [
            ("ORF1ab", 266, 21555),
            ("S", 21563, 25384),
            ("ORF3a", 25393, 26220),
            ("E", 26245, 26472),
            ("M", 26523, 27191),
            ("ORF6", 27202, 27387),
            ("ORF7a", 27394, 27759),
            ("ORF7b", 27756, 27887),
            ("ORF8", 27894, 28259),
            ("N", 28274, 29533),
            ("ORF10", 29558, 29674),
        ]
    elif pathogen.lower() == "rsva":
        gene_annotations = [
            ("NS1", 1, 389),
            ("NS2", 390, 783),
            ("N", 784, 1605),
            ("P", 1606, 2526),
            ("M", 2527, 3378),
            ("SH", 3379, 3723),
            ("G", 3724, 5383),
            ("F", 5384, 7130),
            ("M2", 7131, 7947),
            ("L", 7948, 15191),
        ]
    elif pathogen.lower() == "rsvb":
        gene_annotations = [
            ("NS1", 1, 389),
            ("NS2", 390, 774),
            ("N", 775, 1572),
            ("P", 1573, 2493),
            ("M", 2494, 3345),
            ("SH", 3346, 3690),
            ("G", 3691, 5330),
            ("F", 5331, 7077),
            ("M2", 7078, 7894),
            ("L", 7895, 15132),
        ]
    else:
        raise ValueError(f"Unknown pathogen: {pathogen}")

    # Colors for genes
    gene_colors = dict(zip(
        [g[0] for g in gene_annotations],
        sns.color_palette("tab10", n_colors=len(gene_annotations))
    ))

    bin_mid_to_label = dict(zip(df["bin_mid"], df["bin_label"]))
    unique_bin_labels = list(df["bin_label"].unique())

    for sample_type in df["type"].unique():
        subset = df[df["type"] == sample_type]

        fig, ax = plt.subplots(figsize=(24, 6))
        sns.boxplot(data=subset, x="bin_label", y="coverage", color="skyblue", fliersize=1, ax=ax)

        median_per_bin = subset.groupby("bin_label")["coverage"].median().reset_index()
        sns.lineplot(data=median_per_bin, x="bin_label", y="coverage", color="red", label="Median", linewidth=2, ax=ax)
        # Identify bins where median coverage < 10
        median_per_bin = subset.groupby("bin_label")["coverage"].median().reset_index()
        sns.lineplot(data=median_per_bin, x="bin_label", y="coverage", color="red", label="Median", linewidth=2, ax=ax)

        low_cov_bins = median_per_bin[median_per_bin["coverage"] < 10]["bin_label"].tolist()

        for bin_label in low_cov_bins:
            x_idx = unique_bin_labels.index(bin_label)
            ax.plot(
                x_idx,
                -100,  # slightly below the y-limit max
                marker='*',
                color='black',
                markersize=10,
                zorder=10
            )

        # Draw gene bars under x-axis
        for gene, start, end in gene_annotations:
            gene_color = gene_colors[gene]

            closest_start = min(bin_mid_to_label.keys(), key=lambda x: abs(x - start))
            closest_end = min(bin_mid_to_label.keys(), key=lambda x: abs(x - end))

            if closest_start in bin_mid_to_label and closest_end in bin_mid_to_label:
                x_start = unique_bin_labels.index(bin_mid_to_label[closest_start])
                x_end = unique_bin_labels.index(bin_mid_to_label[closest_end])

                rect = Rectangle(
                    (x_start, -800),
                    width=(x_end - x_start + 1),
                    height=300,
                    facecolor=gene_color,
                    edgecolor='black',
                    linewidth=0.5
                )
                ax.add_patch(rect)

                ax.text(
                    x=(x_start + x_end) / 2,
                    y=-650,
                    s=gene,
                    ha='center',
                    va='center',
                    fontsize=7,
                    rotation=0,
                    color='black'
                )

        ax.set_xticks(np.linspace(0, len(unique_bin_labels) - 1, num=20, dtype=int))
        ax.set_xticklabels([unique_bin_labels[i] for i in ax.get_xticks()], rotation=90, fontsize=6)
        ax.set_xlabel("Genome Position (bin midpoints)")
        ax.set_ylabel("Coverage")
        ax.set_title(f"{pathogen.upper()} Coverage - {sample_type.upper()} Samples\nBoxplots with Median & Gene Annotations")
        ax.set_ylim(-1000, 3000)
        ax.set_yticks(np.arange(300, 3000, 300))

        fig.subplots_adjust(bottom=0.3)

        gene_legend_handles = [Patch(color=gene_colors[gene], label=gene) for gene in gene_colors]
        ax.legend(handles=gene_legend_handles, title="Genes", loc="upper right", fontsize='small', title_fontsize='small')

        output_file = f"{pathogen}_coverage_boxplot_{sample_type}.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot to: {output_file}")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot genome coverage boxplots with gene annotations.")
    parser.add_argument('--input', required=True, help="Path to TSV file with coverage data")
    parser.add_argument('--pathogen', required=True, choices=["sars-cov-2", "rsva", "rsvb"], help="Pathogen type")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    plot_coverage_with_annotations(args.input, pathogen=args.pathogen)

if __name__ == "__main__":
    main()