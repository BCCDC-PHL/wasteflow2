#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from matplotlib.patches import Rectangle
import argparse
import sys
import re

# Import the Metadata parser
from metadata_parser import Metadata

# ----------------------
# Load functions
# ----------------------
def load_mosdepth(file_path):
    df = pd.read_csv(file_path, sep="\t", header=0)
    df_wide = df.pivot_table(index=["chrom", "start", "end"],
                             columns="sample",
                             values="coverage")
    df_wide = df_wide.reset_index()
    cov_cols = df_wide.columns[3:]
    cov_cols_clean = [c.split('_')[0] for c in cov_cols]
    df_wide.rename(columns=dict(zip(cov_cols, cov_cols_clean)), inplace=True)
    return df_wide


def load_bed(file_path):
    bed = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    return bed


# ----------------------
# Normalize coverage sample names
# ----------------------
def normalize_sample_name(s):
    match = re.match(r"(WW\d{2}-\d{3,4})", s)
    return match.group(1) if match else s


# ----------------------
# Compute coverage
# ----------------------
def compute_gene_coverage(mosdepth_df, bed_df, sample_order):
    # Normalize sample names in mosdepth_df columns
    samples_in_data = mosdepth_df.columns[3:]
    normalized_cols = [normalize_sample_name(c) for c in samples_in_data]

    mosdepth_df = mosdepth_df.copy()
    mosdepth_df.columns = list(mosdepth_df.columns[:3]) + normalized_cols

    # Deduplicate columns if normalization caused overlap (keep first)
    mosdepth_df = mosdepth_df.loc[:, ~mosdepth_df.columns.duplicated()]

    gene_cov = {}
    gene_breadth = {}
    
    for _, row in bed_df.iterrows():
        cov_rows = mosdepth_df[(mosdepth_df['start'] < row['end']) & (mosdepth_df['end'] > row['start'])]
        if not cov_rows.empty:
            cov_mean = cov_rows[normalized_cols].mean()
            cov_mean = cov_mean.reindex(sample_order, fill_value=0)
            gene_cov[row['gene']] = cov_mean
            breadth = (cov_rows[normalized_cols] >= 5).sum() / len(cov_rows)
            breadth = breadth.reindex(sample_order, fill_value=0)
            gene_breadth[row['gene']] = breadth
        else:
            gene_cov[row['gene']] = pd.Series([0]*len(sample_order), index=sample_order)
            gene_breadth[row['gene']] = pd.Series([0]*len(sample_order), index=sample_order)
    
    df_cov = pd.DataFrame.from_dict(gene_cov, orient='index').sort_index()
    df_breadth = pd.DataFrame.from_dict(gene_breadth, orient='index').sort_index()
    return df_cov, df_breadth


def make_empty_gene_df(bed_file, sample_list):
    """Make empty coverage & breadth DataFrames if mosdepth not provided."""
    bed = load_bed(bed_file)
    df_cov = pd.DataFrame(0, index=bed['gene'].unique(), columns=sample_list)
    df_breadth = pd.DataFrame(0, index=bed['gene'].unique(), columns=sample_list)
    return df_cov, df_breadth


# ----------------------
# Plot heatmap
# ----------------------
def plot_heatmap_with_breadth(df_cov, df_breadth, pathogen_name, ax=None, vmin=5, vmax=100, bold_threshold=0.75, x_labels=None):
    df_plot = df_cov.clip(lower=vmin, upper=vmax)
    cmap = LinearSegmentedColormap.from_list("white_seagreen", ["white", "lightgreen", "seagreen"])
    
    sns.heatmap(df_plot, cmap=cmap, linewidths=0.5, linecolor='gray',
                ax=ax, vmin=vmin, vmax=vmax, annot=False, cbar_kws={'label':'Mean coverage (X)'})
    
    ax.set_title(f"{pathogen_name}")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    if x_labels is not None:
        ax.set_xticklabels(x_labels, rotation=90, ha='center')
    else:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    nrows, ncols = df_plot.shape
    for i in range(nrows):
        for j in range(ncols):
            if df_breadth.iloc[i,j] >= bold_threshold:
                rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=2.5)
                ax.add_patch(rect)


# ----------------------
# Main
# ----------------------
def main():
    parser = argparse.ArgumentParser(description="Multi-panel gene-level coverage heatmap with metadata filtering")
    parser.add_argument("--sars", required=False, help="SARS-CoV-2 coverage file (mosdepth)")
    parser.add_argument("--rsvA", required=False, help="RSV-A coverage file (mosdepth)")
    parser.add_argument("--rsvB", required=False, help="RSV-B coverage file (mosdepth)")
    parser.add_argument("--sars_bed", required=True, help="SARS-CoV-2 gene BED file")
    parser.add_argument("--rsvA_bed", required=True, help="RSV-A gene BED file")
    parser.add_argument("--rsvB_bed", required=True, help="RSV-B gene BED file")
    parser.add_argument("--metadata", required=False, help="Optional metadata file (tsv)")
    parser.add_argument("--filter", required=False, nargs='+',
                        help="Optional filters to apply: e.g., season=2025_2026 period=in_season extraction_source=research")
    parser.add_argument("--out", required=True, help="Output figure file name (png)")
    args = parser.parse_args()
    
    # ----------------------
    # Determine sample order using Metadata parser
    # ----------------------
    if args.metadata:
        md = Metadata(args.metadata)
        if args.filter:
            filters = dict(f.split("=") for f in args.filter)
            meta_filtered = md.filter(**filters)
        else:
            meta_filtered = md.df.copy()

        if meta_filtered.empty:
            print("⚠️ No samples found after applying filters.")
            sys.exit(0)

        # Ensure collection_date is datetime
        meta_filtered['collection_date'] = pd.to_datetime(meta_filtered['collection_date'])

        # Define custom site order for sorting
        site_order = ["NHA", "IHA", "MV", "VIHA"]
        meta_filtered['site_prefix'] = meta_filtered['collection_site'].str.split('-').str[0]
        meta_filtered['site_prefix'] = pd.Categorical(meta_filtered['site_prefix'], categories=site_order, ordered=True)

        # Sort first by date, then by site order
        meta_sorted = meta_filtered.sort_values(['collection_date', 'site_prefix', 'collection_site'])

        sample_list = meta_sorted["sample_id"].tolist()
        x_labels = meta_sorted.apply(
            lambda r: f"{r['collection_site']}\n{r['collection_date'].strftime('%m/%d')}", axis=1
        ).tolist()
    else:
        # Default: infer from any coverage file
        if args.sars:
            sample_list = load_mosdepth(args.sars).columns[3:].tolist()
        elif args.rsvA:
            sample_list = load_mosdepth(args.rsvA).columns[3:].tolist()
        elif args.rsvB:
            sample_list = load_mosdepth(args.rsvB).columns[3:].tolist()
        else:
            raise ValueError("No metadata or coverage files provided → cannot infer sample list")
        x_labels = sample_list

    # ----------------------
    # SARS-CoV-2
    # ----------------------
    if args.sars:
        df_sars_cov, df_sars_breadth = compute_gene_coverage(load_mosdepth(args.sars), load_bed(args.sars_bed), sample_list)
    else:
        df_sars_cov, df_sars_breadth = make_empty_gene_df(args.sars_bed, sample_list)

    # RSV-A
    if args.rsvA:
        df_rsvA_cov, df_rsvA_breadth = compute_gene_coverage(load_mosdepth(args.rsvA), load_bed(args.rsvA_bed), sample_list)
    else:
        df_rsvA_cov, df_rsvA_breadth = make_empty_gene_df(args.rsvA_bed, sample_list)

    # RSV-B
    if args.rsvB:
        df_rsvB_cov, df_rsvB_breadth = compute_gene_coverage(load_mosdepth(args.rsvB), load_bed(args.rsvB_bed), sample_list)
    else:
        df_rsvB_cov, df_rsvB_breadth = make_empty_gene_df(args.rsvB_bed, sample_list)

    # ----------------------
    # Plot multi-panel figure
    # ----------------------
    fig, axes = plt.subplots(3, 1, figsize=(max(15, len(sample_list)*0.5), 12))
    plot_heatmap_with_breadth(df_sars_cov, df_sars_breadth, "SARS-CoV-2", ax=axes[0], x_labels=x_labels)
    plot_heatmap_with_breadth(df_rsvA_cov, df_rsvA_breadth, "RSV-A", ax=axes[1], x_labels=x_labels)
    plot_heatmap_with_breadth(df_rsvB_cov, df_rsvB_breadth, "RSV-B", ax=axes[2], x_labels=x_labels)

    fig.text(0.5, 0.02, "Cells with black borders indicate genes where ≥75% of the gene has ≥5× coverage",
             ha='center', fontsize=10)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(args.out, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()