#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot QC results for viral and control oligos"
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


def main():
    args = parse_args()

    # Read the concatenated TSV file
    df = pd.read_csv(
        args.input_file,
        sep="\t",
        header=None,
        names=["sample", "ref", "length", "mapped1", "mapped2"],
    )

    # Sum mapped reads
    df["total_reads"] = df["mapped1"] + df["mapped2"]

    # Identify viral refs = everything except control and other
    viral_refs = sorted(set(df["ref"]) - {args.control_ref, "*"})

    # Pivot to wide format
    reads = df.pivot(index="sample", columns="ref", values="total_reads").fillna(0)

    # Calculate total reads per sample
    reads["total"] = reads.sum(axis=1)

    # Calculate proportions
    for ref in viral_refs:
        reads[f"{ref}_prop"] = reads[ref] / reads["total"]
    reads["Control_prop"] = reads[args.control_ref] / reads["total"]
    reads["Unmapped_prop"] = reads["*"] / reads["total"]

    # Build proportion dataframe for plotting
    plot_df = reads[
        [f"{ref}_prop" for ref in viral_refs] + ["Control_prop", "Unmapped_prop"]
    ]

    # Plot stacked bar
    ax = plot_df.plot(kind="bar", stacked=True, figsize=(14, 6), colormap="tab20")

    # Add total reads on top of each bar
    for i, total in enumerate(reads["total"]):
        ax.text(i, 1.05, f"{int(total):,}", ha="center", va="bottom", fontsize=8)

    plt.ylabel("Proportion of reads")
    plt.title(
        "Proportion of Viral, Control, and Unmapped Reads per Sample\n(Total reads shown above each bar)"
    )
    plt.xticks(rotation=90, ha="right")
    # plt.legend([*viral_refs, "Control", "Unmapped"])
    plt.legend(
        [*viral_refs, "Control", "Unmapped"], loc="center left", bbox_to_anchor=(1, 0.5)
    )
    plt.tight_layout()

    # Save plot
    out_plot = f"{args.out_prefix}_stacked_bar.png"
    plt.savefig(out_plot, dpi=300)
    plt.close()

    # Save table of absolute and proportion reads
    out_table = f"{args.out_prefix}_read_counts.tsv"
    reads.to_csv(out_table, sep="\t")

    print(f"[INFO] Found viral references: {', '.join(viral_refs)}")
    print(f"[INFO] Plot saved to {out_plot}")
    print(f"[INFO] Table saved to {out_table}")


if __name__ == "__main__":
    main()
