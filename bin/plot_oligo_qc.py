#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

# Regex to detect hg38 human contigs
HUMAN_CHR_REGEX = re.compile(r"^chr([0-9]+|X|Y|M|Un|KI|GL).*")

# List of all Influenza accessions in your panel
INFLUENZA_A_ACCESSIONS = [
    "NC_026438.1", "NC_026435.1", "NC_026437.1", "NC_026433.1",
    "NC_026436.1", "NC_026434.1", "NC_026431.1", "NC_026432.1",
    "NC_007373.1", "NC_007372.1", "NC_007371.1", "CY163680.1",
    "NC_007369.1", "CY114383.1", "NC_007367.1", "NC_007370.1",
    "NC_007357.1", "NC_007358.1", "NC_007359.1", "NC_007362.1",
    "NC_007360.1", "NC_007361.1", "NC_007363.1", "NC_007364.1"
]
INFLUENZA_B_VIC_ACCESSIONS = [
    "CY115158.1", "CY115157.1", "CY115156.1", "KX058884.1", 
    "CY115154.1", "CY073894.1", "CY115152.1", "CY115155.1" 
]

SARS_COV_ACCESSION = "MN908947.3"
RSV_A_ACCESSION = "PP109421.1"
RSV_B_ACCESSION = "OP975389.1"
MEASLES_ACCESSION = "NC_001498.1"



def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot QC results for viral panel including human and influenza"
    )
    parser.add_argument(
        "--input_file", required=True, help="Input TSV file with concatenated idxstats"
    )
    parser.add_argument(
        "--control_ref", default="CTRL_OLIGO", help="Control oligo reference ID"
    )
    parser.add_argument(
        "--out_prefix", default="qc_results", help="Prefix for output files"
    )
    return parser.parse_args()


def classify_ref(ref, control_ref):
    if ref == "*":
        return "Unmapped"
    if ref == control_ref:
        return "Control"
    if ref in INFLUENZA_A_ACCESSIONS:
        return "Influenza A"
    if ref in INFLUENZA_B_VIC_ACCESSIONS:
        return "Influenza B (Victoria)"
    if HUMAN_CHR_REGEX.match(ref):
        return "Human"
    if ref == SARS_COV_ACCESSION:
        return "SARS-CoV-2"
    if ref == RSV_A_ACCESSION:
        return "RSV-A"
    if ref == RSV_B_ACCESSION:
        return "RSV-B"
    if ref == MEASLES_ACCESSION:
        return "Measles"
    return ref


def main():
    args = parse_args()

    df = pd.read_csv(
        args.input_file,
        sep="\t",
        header=None,
        names=["sample", "ref", "length", "mapped1", "mapped2"],
    )

    df["total_reads"] = df["mapped1"] + df["mapped2"]
    df["category"] = df["ref"].apply(lambda x: classify_ref(x, args.control_ref))

    reads = (
        df.groupby(["sample", "category"])["total_reads"]
        .sum()
        .unstack(fill_value=0)
    )

    # Remove ".bam" suffix from sample names
    reads.index = reads.index.str.replace(r"\.bam$", "", regex=True)

    reads["total"] = reads.sum(axis=1)

    for col in reads.columns:
        if col != "total":
            reads[f"{col}_prop"] = reads[col] / reads["total"]

    category_order = [
        "SARS-CoV-2",
        "RSV-B",
        "RSV-A",
        "Influenza A",
        "Influenza B (Victoria)",
        "Measles",
        "Control",
        "Human",
        "Unmapped",
    ]
    category_order = [c for c in category_order if c in reads.columns]

    plot_df = reads[[f"{c}_prop" for c in category_order]]

    colors = {
        "SARS-CoV-2": "#1f77b4",
        "RSV-B": "#ff7f0e",
        "RSV-A": "#2ca02c",
        "Influenza A": "#d62728",
        "Influenza B (Victoria)": "#9467bd",
        "Measles": "#8c564b",
        "Control": "#e377c2",
        "Human": "#bcbd22",
        "Unmapped": "#7f7f7f",
    }
    plot_colors = [colors[c] for c in category_order]

    fig, ax = plt.subplots(figsize=(14, 6))
    plot_df.plot(kind="bar", stacked=True, ax=ax, color=plot_colors)

    # Total read labels (unchanged)
    for i, total in enumerate(reads["total"]):
        ax.text(
            i,
            1.05,
            f"{int(total):,}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_ylabel("Proportion of reads")
    plt.xticks(rotation=90, ha="right")
    plt.yticks(rotation=90, ha="right")
    ax.legend(category_order, loc="center left", bbox_to_anchor=(1, 0.5))

    # ---- Bottom title ----
    #fig.text(
    #    0.5,
    #   0.01,
    #    "Proportion of Viral, Influenza, Control, Human, and Unmapped Reads per Sample\n"
    #    "(Total reads shown above each bar)",
    #    ha="center",
    #    va="bottom",
    #    fontsize=12,
    #)

    # Make space for bottom title
    #plt.subplots_adjust(bottom=0.22)

    out_plot = f"{args.out_prefix}_stacked_bar.png"
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    plt.close()

    out_table = f"{args.out_prefix}_read_counts.tsv"
    reads.to_csv(out_table, sep="\t")

    print(f"[INFO] Plot saved to {out_plot}")
    print(f"[INFO] Table saved to {out_table}")


if __name__ == "__main__":
    main()