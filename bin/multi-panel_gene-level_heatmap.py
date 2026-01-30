#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from matplotlib.patches import Rectangle
import argparse
import sys
import re
import os

from metadata_parser import Metadata

# ----------------------
# Load functions
# ----------------------
def load_mosdepth(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_path} does not exist.")
    df = pd.read_csv(file_path, sep="\t", header=0)
    if df.empty:
        raise ValueError(f"{file_path} is empty.")
    if "coverage" not in df.columns:
        raise ValueError(f"{file_path} does not contain 'coverage'. Columns: {df.columns.tolist()}")
    return df

def load_bed(file_path):
    return pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "gene"],
        na_filter=False
    )

def normalize_sample_name(s, influenza=False):
    """Extract base sample name without segment/subtype suffixes"""
    match = re.match(r"(WW\d{2}-\d{3,4})", s)
    if match:
        return match.group(1)
    match2 = re.match(r"(Undetermined_S\d+)", s)
    if match2:
        return match2.group(1)
    return s

def extract_segment_name(sample_col, influenza=False):
    """Extract segment name from sample column for influenza"""
    if not influenza:
        return None
    # For influenza samples like: WW25-1508_S5_H5N1_M
    # Extract the segment part (M, PB1, PB2, etc.)
    parts = sample_col.split('_')
    if len(parts) >= 3:
        return parts[-1]  # Last part is segment
    return None

def compute_gene_coverage_influenza(mosdepth_df, bed_df, sample_order=None):
    """
    Special handling for influenza where each sample has multiple segment columns
    """
    mosdepth_df = mosdepth_df.copy()
    
    # Create a mapping: chrom -> segment_name
    chrom_to_segment = dict(zip(bed_df['chrom'], bed_df['gene']))
    
    # Add segment column to mosdepth data based on chrom
    mosdepth_df['segment'] = mosdepth_df['chrom'].map(chrom_to_segment)
    
    # Extract base sample name and create sample column
    mosdepth_df['base_sample'] = mosdepth_df['sample'].apply(lambda x: normalize_sample_name(x, influenza=True))
    
    # Group by segment, base_sample, and genomic regions, then calculate mean coverage
    # This handles multiple entries per segment per sample
    grouped = mosdepth_df.groupby(['segment', 'base_sample', 'chrom', 'start', 'end'])['coverage'].mean().reset_index()
    
    # Now compute per-segment, per-sample statistics
    gene_cov = {}
    gene_breadth = {}
    
    for segment in bed_df['gene'].unique():
        segment_data = grouped[grouped['segment'] == segment]
        
        if segment_data.empty:
            # No data for this segment
            if sample_order:
                gene_cov[segment] = pd.Series(0, index=sample_order)
                gene_breadth[segment] = pd.Series(0, index=sample_order)
            continue
        
        # Pivot to get samples as columns
        pivot = segment_data.pivot_table(
            index=['start', 'end'], 
            columns='base_sample', 
            values='coverage',
            aggfunc='mean'
        )
        
        # Ensure all samples in sample_order are present
        if sample_order:
            for s in sample_order:
                if s not in pivot.columns:
                    pivot[s] = 0
            pivot = pivot[sample_order]
        
        # Calculate mean coverage and breadth
        mean_cov = pivot.mean()
        breadth = (pivot >= 5).sum() / len(pivot) if len(pivot) > 0 else 0
        
        gene_cov[segment] = mean_cov
        gene_breadth[segment] = breadth
    
    # Create dataframes with segment order from BED
    gene_order = bed_df['gene'].tolist()
    df_cov = pd.DataFrame.from_dict(gene_cov, orient='index')
    df_breadth = pd.DataFrame.from_dict(gene_breadth, orient='index')
    
    # Reindex to match BED file order
    df_cov = df_cov.reindex(gene_order)
    df_breadth = df_breadth.reindex(gene_order)
    
    # Fill NaN with 0
    df_cov = df_cov.fillna(0)
    df_breadth = df_breadth.fillna(0)
    
    return df_cov, df_breadth

def compute_gene_coverage(mosdepth_df, bed_df, sample_order=None, influenza=False):
    """Standard coverage computation for single-segment viruses (SARS-CoV-2, RSV)"""
    if influenza:
        return compute_gene_coverage_influenza(mosdepth_df, bed_df, sample_order)
    
    # Original logic for non-influenza
    df_wide = mosdepth_df.pivot_table(
        index=["chrom", "start", "end"], columns="sample", values="coverage"
    ).reset_index()
    
    cov_cols = df_wide.columns[3:]
    cov_cols_clean = [normalize_sample_name(c) for c in cov_cols]
    df_wide.rename(columns=dict(zip(cov_cols, cov_cols_clean)), inplace=True)
    
    # Handle duplicate columns by averaging
    df_numeric = df_wide.iloc[:, 3:]
    if df_numeric.columns.duplicated().any():
        df_numeric = df_numeric.T.groupby(level=0).mean().T
    
    # Ensure all samples present
    if sample_order:
        for s in sample_order:
            if s not in df_numeric.columns:
                df_numeric[s] = 0
        df_numeric = df_numeric[sample_order]
    
    df_wide_clean = pd.concat([df_wide.iloc[:, :3], df_numeric], axis=1)
    
    gene_cov = {}
    gene_breadth = {}
    for _, row in bed_df.iterrows():
        cov_rows = df_wide_clean[
            (df_wide_clean["start"] < row["end"]) & 
            (df_wide_clean["end"] > row["start"])
        ]
        if not cov_rows.empty:
            cov_mean = cov_rows.iloc[:, 3:].mean()
            breadth = (cov_rows.iloc[:, 3:] >= 5).sum() / len(cov_rows)
        else:
            cov_mean = pd.Series(0, index=df_numeric.columns)
            breadth = pd.Series(0, index=df_numeric.columns)
        gene_cov[row["gene"]] = cov_mean
        gene_breadth[row["gene"]] = breadth
    
    gene_order = bed_df["gene"].tolist()
    df_cov = pd.DataFrame.from_dict(gene_cov, orient="index").reindex(gene_order)
    df_breadth = pd.DataFrame.from_dict(gene_breadth, orient="index").reindex(gene_order)
    
    return df_cov, df_breadth

def make_empty_gene_df(bed_file, sample_list):
    bed = load_bed(bed_file)
    gene_order = bed["gene"].tolist()
    df_cov = pd.DataFrame(0, index=gene_order, columns=sample_list)
    df_breadth = pd.DataFrame(0, index=gene_order, columns=sample_list)
    return df_cov, df_breadth

def generate_x_labels(df_cov, metadata_df=None):
    """Generate x-axis labels matching the columns of df_cov"""
    cols = df_cov.columns.tolist()
    if metadata_df is not None:
        label_map = {
            row["sample_id"]: f"{row['collection_site']}\n{pd.to_datetime(row['collection_date']).strftime('%m/%d')}"
            for _, row in metadata_df.iterrows()
        }
        labels = [label_map.get(c, c) for c in cols]
    else:
        labels = cols
    return labels

def plot_heatmap_with_breadth(df_cov, df_breadth, pathogen_name, ax=None, vmin=5, vmax=100, bold_threshold=0.75, x_labels=None):
    df_plot = df_cov.clip(lower=vmin, upper=vmax)
    cmap = LinearSegmentedColormap.from_list("white_seagreen", ["white", "lightgreen", "seagreen"])
    sns.heatmap(df_plot, cmap=cmap, linewidths=0.5, linecolor="gray", ax=ax,
                vmin=vmin, vmax=vmax, annot=False, cbar_kws={"label": "Mean coverage (X)"})
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
    ax.set_title(f"{pathogen_name}")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes" if "RSV" in pathogen_name or "SARS" in pathogen_name else "Segments")
    if x_labels is not None:
        if len(x_labels) == df_plot.shape[1]:
            ax.set_xticklabels(x_labels, rotation=90, ha="center")
        else:
            ax.set_xticklabels(df_plot.columns, rotation=90, ha="center")
    
    nrows, ncols = df_plot.shape
    for i in range(nrows):
        for j in range(ncols):
            if df_breadth.iloc[i, j] >= bold_threshold:
                rect = Rectangle((j, i), 1, 1, fill=False, edgecolor="black", lw=2.5)
                ax.add_patch(rect)

def load_or_empty(cov_file, bed_file, sample_list, influenza=False):
    if cov_file and os.path.exists(cov_file):
        mosdepth_df = load_mosdepth(cov_file)
        bed_df = load_bed(bed_file)
        return compute_gene_coverage(mosdepth_df, bed_df, sample_list, influenza=influenza)
    return make_empty_gene_df(bed_file, sample_list)

def write_segment_coverage_tsv(df_breadth, output_file, level_name="Segment"):
    df = df_breadth.copy()
    if isinstance(df.index, pd.MultiIndex):
        df_long = df.reset_index().melt(
            id_vars=list(df.index.names),
            var_name="Sample",
            value_name="percent_covered"
        )
    else:
        df_long = df.reset_index().melt(
            id_vars="index",
            var_name="Sample",
            value_name="percent_covered"
        ).rename(columns={"index": level_name})
    df_long["percent_covered"] = (df_long["percent_covered"] * 100).round(1)
    df_long.to_csv(output_file, sep="\t", index=False)

# ----------------------
# Main
# ----------------------
def main():
    parser = argparse.ArgumentParser(description="Unified multi-panel coverage heatmap")
    parser.add_argument("--sars", help="SARS-CoV-2 coverage file (mosdepth)")
    parser.add_argument("--rsvA", help="RSV-A coverage file (mosdepth)")
    parser.add_argument("--rsvB", help="RSV-B coverage file (mosdepth)")
    parser.add_argument("--h1n1", help="H1N1 coverage file (mosdepth)")
    parser.add_argument("--h3n2", help="H3N2 coverage file (mosdepth)")
    parser.add_argument("--h5n1", help="H5N1 coverage file (mosdepth)")
    parser.add_argument("--sars_bed", required=True)
    parser.add_argument("--rsvA_bed", required=True)
    parser.add_argument("--rsvB_bed", required=True)
    parser.add_argument("--h1n1_bed", required=True)
    parser.add_argument("--h3n2_bed", required=True)
    parser.add_argument("--h5n1_bed", required=True)
    parser.add_argument("--metadata", help="Optional metadata file (tsv)")
    parser.add_argument("--filter", nargs="+")
    parser.add_argument("--out_sarsrsv", required=True, help="Output figure for SARS/RSV")
    parser.add_argument("--out_influenza", required=True, help="Output figure for Influenza")
    parser.add_argument("--out_sarsrsv_tsv", help="Output TSV for SARS/RSV coverage")
    parser.add_argument("--out_influenza_tsv", help="Output TSV for Influenza coverage")
    args = parser.parse_args()

    # ----------------------
    # Determine sample order
    # ----------------------
    if args.metadata:
        md = Metadata(args.metadata)
        meta_filtered = md.df.copy()
        if args.filter:
            filters = dict(f.split("=") for f in args.filter)
            meta_filtered = md.filter(**filters)
        if meta_filtered.empty:
            print("⚠️ No samples found after applying filters.")
            sys.exit(0)
        meta_filtered["collection_date"] = pd.to_datetime(meta_filtered["collection_date"])
        site_order = ["NHA", "IHA", "MV", "VIHA"]
        meta_filtered["site_prefix"] = meta_filtered["collection_site"].str.split("-").str[0]
        meta_filtered["site_prefix"] = pd.Categorical(meta_filtered["site_prefix"], categories=site_order, ordered=True)
        meta_sorted = meta_filtered.sort_values(["collection_date", "site_prefix", "collection_site"])
        sample_list = meta_sorted["sample_id"].tolist()
    else:
        # Extract sample list from first available file
        sample_list = None
        for cov in [args.sars, args.rsvA, args.rsvB, args.h1n1, args.h3n2, args.h5n1]:
            if cov and os.path.exists(cov):
                df_temp = load_mosdepth(cov)
                # For influenza, extract base sample names
                is_influenza = any(x in cov for x in ['h1n1', 'h3n2', 'h5n1'])
                samples = df_temp['sample'].unique()
                sample_list = sorted(list(set([normalize_sample_name(s, influenza=is_influenza) for s in samples])))
                break
        meta_sorted = None

    # ----------------------
    # Load SARS/RSV
    # ----------------------
    df_sars_cov, df_sars_breadth = load_or_empty(args.sars, args.sars_bed, sample_list, influenza=False)
    df_rsvA_cov, df_rsvA_breadth = load_or_empty(args.rsvA, args.rsvA_bed, sample_list, influenza=False)
    df_rsvB_cov, df_rsvB_breadth = load_or_empty(args.rsvB, args.rsvB_bed, sample_list, influenza=False)

    x_labels_sars = generate_x_labels(df_sars_cov, meta_sorted)
    x_labels_rsvA = generate_x_labels(df_rsvA_cov, meta_sorted)
    x_labels_rsvB = generate_x_labels(df_rsvB_cov, meta_sorted)

    if args.out_sarsrsv_tsv:
        df_sarsrsv_breadth = pd.concat(
            [df_sars_breadth, df_rsvA_breadth, df_rsvB_breadth],
            keys=["SARS-CoV-2","RSV-A","RSV-B"],
            names=["Pathogen","Gene"]
        )
        write_segment_coverage_tsv(df_sarsrsv_breadth, args.out_sarsrsv_tsv, level_name="Gene")

    fig, axes = plt.subplots(3,1, figsize=(max(15, max(df_sars_cov.shape[1], df_rsvA_cov.shape[1], df_rsvB_cov.shape[1])*0.5), 12))
    plot_heatmap_with_breadth(df_sars_cov, df_sars_breadth, "SARS-CoV-2", ax=axes[0], x_labels=x_labels_sars)
    plot_heatmap_with_breadth(df_rsvA_cov, df_rsvA_breadth, "RSV-A", ax=axes[1], x_labels=x_labels_rsvA)
    plot_heatmap_with_breadth(df_rsvB_cov, df_rsvB_breadth, "RSV-B", ax=axes[2], x_labels=x_labels_rsvB)
    fig.text(0.5,0.02,"Cells with black borders indicate genes where ≥75% of the gene has ≥5× coverage", ha="center", fontsize=10)
    plt.tight_layout(rect=[0,0.05,1,1])
    plt.savefig(args.out_sarsrsv, dpi=300)
    plt.close(fig)

    # ----------------------
    # Load Influenza
    # ----------------------
    df_h1n1_cov, df_h1n1_breadth = load_or_empty(args.h1n1, args.h1n1_bed, sample_list, influenza=True)
    df_h3n2_cov, df_h3n2_breadth = load_or_empty(args.h3n2, args.h3n2_bed, sample_list, influenza=True)
    df_h5n1_cov, df_h5n1_breadth = load_or_empty(args.h5n1, args.h5n1_bed, sample_list, influenza=True)

    x_labels_h1n1 = generate_x_labels(df_h1n1_cov, meta_sorted)
    x_labels_h3n2 = generate_x_labels(df_h3n2_cov, meta_sorted)
    x_labels_h5n1 = generate_x_labels(df_h5n1_cov, meta_sorted)

    if args.out_influenza_tsv:
        df_influenza_breadth = pd.concat(
            [df_h1n1_breadth, df_h3n2_breadth, df_h5n1_breadth],
            keys=["H1N1","H3N2","H5N1"],
            names=["Pathogen","Segment"]
        )
        write_segment_coverage_tsv(df_influenza_breadth, args.out_influenza_tsv, level_name="Segment")

    fig, axes = plt.subplots(3,1, figsize=(max(15, max(df_h1n1_cov.shape[1], df_h3n2_cov.shape[1], df_h5n1_cov.shape[1])*0.5), 12))
    plot_heatmap_with_breadth(df_h1n1_cov, df_h1n1_breadth, "Influenza A (H1N1)", ax=axes[0], bold_threshold=0.75, x_labels=x_labels_h1n1)
    plot_heatmap_with_breadth(df_h3n2_cov, df_h3n2_breadth, "Influenza A (H3N2)", ax=axes[1], bold_threshold=0.75, x_labels=x_labels_h3n2)
    plot_heatmap_with_breadth(df_h5n1_cov, df_h5n1_breadth, "Influenza A (H5N1)", ax=axes[2], bold_threshold=0.75, x_labels=x_labels_h5n1)
    fig.text(0.5,0.02,"Cells with black borders indicate segments where ≥75% of the segment has ≥5× coverage", ha="center", fontsize=10)
    plt.tight_layout(rect=[0,0.05,1,1])
    plt.savefig(args.out_influenza, dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    main()