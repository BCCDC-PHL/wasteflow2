#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from matplotlib.patches import Rectangle
import argparse

# ----------------------
# Load functions
# ----------------------
def load_mosdepth(file_path):
    df = pd.read_csv(file_path, sep="\t", header=0)
    df_wide = df.pivot_table(index=["chrom", "start", "end"],
                             columns="sample",
                             values="coverage")
    df_wide = df_wide.reset_index()
    # Normalize coverage columns: remove "_SXX" suffix
    cov_cols = df_wide.columns[3:]
    cov_cols_clean = [c.split('_')[0] for c in cov_cols]
    df_wide.rename(columns=dict(zip(cov_cols, cov_cols_clean)), inplace=True)
    return df_wide

def load_bed(file_path):
    bed = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    return bed

def compute_gene_coverage(mosdepth_df, bed_df, sample_order):
    samples_in_data = mosdepth_df.columns[3:]
    gene_cov = {}
    gene_breadth = {}
    
    for _, row in bed_df.iterrows():
        cov_rows = mosdepth_df[(mosdepth_df['start'] < row['end']) & (mosdepth_df['end'] > row['start'])]
        if not cov_rows.empty:
            cov_mean = cov_rows[samples_in_data].mean()
            cov_mean = cov_mean.reindex(sample_order, fill_value=0)
            gene_cov[row['gene']] = cov_mean
            
            # Fraction of gene ≥5x
            breadth = (cov_rows[samples_in_data] >= 5).sum() / len(cov_rows)
            breadth = breadth.reindex(sample_order, fill_value=0)
            gene_breadth[row['gene']] = breadth
        else:
            gene_cov[row['gene']] = pd.Series([0]*len(sample_order), index=sample_order)
            gene_breadth[row['gene']] = pd.Series([0]*len(sample_order), index=sample_order)
    
    df_cov = pd.DataFrame.from_dict(gene_cov, orient='index').sort_index()
    df_breadth = pd.DataFrame.from_dict(gene_breadth, orient='index').sort_index()
    return df_cov, df_breadth

# ----------------------
# Plot heatmap
# ----------------------
def plot_heatmap_with_breadth(df_cov, df_breadth, pathogen_name, ax=None, vmin=5, vmax=100, bold_threshold=0.8, x_labels=None):
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
    # Bold rectangles for cells meeting threshold
    for i in range(nrows):
        for j in range(ncols):
            if df_breadth.iloc[i,j] >= bold_threshold:
                rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=2.5)
                ax.add_patch(rect)

# ----------------------
# Main
# ----------------------
def main():
    parser = argparse.ArgumentParser(description="Multi-panel gene-level coverage heatmap with optional metadata sorting")
    parser.add_argument("--sars", required=True, help="SARS-CoV-2 coverage file (mosdepth)")
    parser.add_argument("--rsvA", required=True, help="RSV-A coverage file (mosdepth)")
    parser.add_argument("--rsvB", required=True, help="RSV-B coverage file (mosdepth)")
    parser.add_argument("--sars_bed", required=True, help="SARS-CoV-2 gene BED file")
    parser.add_argument("--rsvA_bed", required=True, help="RSV-A gene BED file")
    parser.add_argument("--rsvB_bed", required=True, help="RSV-B gene BED file")
    parser.add_argument("--metadata", required=False, help="Optional metadata file (tsv) with Sample_ID, Batch_name, Patient First Name, Collection Date")
    parser.add_argument("--batch", required=False, help="Batch name to filter samples in metadata")
    parser.add_argument("--out", required=True, help="Output figure file name (png)")
    args = parser.parse_args()
    
    # ----------------------
    # Determine sample order
    # ----------------------
    if args.metadata and args.batch:
        meta = pd.read_csv(args.metadata, sep="\t")
        meta = meta[meta['Batch_name']==args.batch].copy()
        meta['Sample_ID_clean'] = meta['Sample_ID'].str.replace(' ', '-')
        meta['Collection Date'] = pd.to_datetime(meta['Collection Date'])
        meta_sorted = meta.sort_values(['Patient First Name', 'Collection Date'])
        sample_list = meta_sorted['Sample_ID_clean'].tolist()
        x_labels = meta_sorted.apply(lambda r: f"{r['Patient First Name']}\n{r['Collection Date'].strftime('%m/%d')}", axis=1).tolist()
    else:
        sample_list = None
        x_labels = None

    # Load coverage files
    df_sars = load_mosdepth(args.sars)
    df_rsvA = load_mosdepth(args.rsvA)
    df_rsvB = load_mosdepth(args.rsvB)

    # If no sample_list provided, just use coverage columns
    if sample_list is None:
        sample_list = df_sars.columns[3:].tolist()
        x_labels = sample_list

    # Compute gene coverage & breadth
    df_sars_cov, df_sars_breadth = compute_gene_coverage(df_sars, load_bed(args.sars_bed), sample_list)
    df_rsvA_cov, df_rsvA_breadth = compute_gene_coverage(df_rsvA, load_bed(args.rsvA_bed), sample_list)
    df_rsvB_cov, df_rsvB_breadth = compute_gene_coverage(df_rsvB, load_bed(args.rsvB_bed), sample_list)

    # Plot multi-panel figure
    fig, axes = plt.subplots(3, 1, figsize=(max(15, len(sample_list)*0.5),12))
    plot_heatmap_with_breadth(df_sars_cov, df_sars_breadth, "SARS-CoV-2", ax=axes[0], x_labels=x_labels)
    plot_heatmap_with_breadth(df_rsvA_cov, df_rsvA_breadth, "RSV-A", ax=axes[1], x_labels=x_labels)
    plot_heatmap_with_breadth(df_rsvB_cov, df_rsvB_breadth, "RSV-B", ax=axes[2], x_labels=x_labels)

    # Add explanatory text
    fig.text(0.5, 0.02, "Cells with black borders indicate genes where ≥80% of the gene has ≥5× coverage",
             ha='center', fontsize=10)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(args.out, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()