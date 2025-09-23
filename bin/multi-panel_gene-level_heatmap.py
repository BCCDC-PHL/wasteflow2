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
    return df_wide

def load_bed(file_path):
    bed = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    return bed

def compute_gene_coverage(mosdepth_df, bed_df, all_samples):
    samples_in_data = mosdepth_df.columns[3:]
    gene_cov = {}
    gene_breadth = {}
    
    for _, row in bed_df.iterrows():
        cov_rows = mosdepth_df[(mosdepth_df['start'] < row['end']) & (mosdepth_df['end'] > row['start'])]
        if not cov_rows.empty:
            cov_mean = cov_rows[samples_in_data].mean()
            cov_mean = cov_mean.reindex(all_samples, fill_value=0)
            gene_cov[row['gene']] = cov_mean
            
            # Fraction of gene ≥5x
            breadth = (cov_rows[samples_in_data] >= 5).sum() / len(cov_rows)
            breadth = breadth.reindex(all_samples, fill_value=0)
            gene_breadth[row['gene']] = breadth
        else:
            gene_cov[row['gene']] = pd.Series([0]*len(all_samples), index=all_samples)
            gene_breadth[row['gene']] = pd.Series([0]*len(all_samples), index=all_samples)
    
    df_cov = pd.DataFrame.from_dict(gene_cov, orient='index').sort_index()
    df_breadth = pd.DataFrame.from_dict(gene_breadth, orient='index').sort_index()
    return df_cov, df_breadth

# ----------------------
# Plot heatmap
# ----------------------
def plot_heatmap_with_breadth(df_cov, df_breadth, pathogen_name, ax=None, vmin=5, vmax=100, bold_threshold=0.8):
    df_plot = df_cov.clip(lower=vmin, upper=vmax)
    cmap = LinearSegmentedColormap.from_list("white_seagreen", ["white", "lightgreen", "seagreen"])
    
    # Plot heatmap
    sns.heatmap(df_plot, cmap=cmap, linewidths=0.5, linecolor='gray',
                ax=ax, vmin=vmin, vmax=vmax, annot=False, cbar_kws={'label':'Mean coverage (X)'})
    
    ax.set_title(f"{pathogen_name}")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    nrows, ncols = df_plot.shape
    
    # Add bold rectangle for cells meeting breadth threshold
    for i in range(nrows):
        for j in range(ncols):
            if df_breadth.iloc[i,j] >= bold_threshold:
                rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=2.5)
                ax.add_patch(rect)

# ----------------------
# Main
# ----------------------
def main():
    parser = argparse.ArgumentParser(description="Multi-panel gene-level coverage heatmap")
    parser.add_argument("--sars", required=True, help="SARS-CoV-2 coverage file (mosdepth)")
    parser.add_argument("--rsvA", required=True, help="RSV-A coverage file (mosdepth)")
    parser.add_argument("--rsvB", required=True, help="RSV-B coverage file (mosdepth)")
    parser.add_argument("--sars_bed", required=True, help="SARS-CoV-2 gene BED file")
    parser.add_argument("--rsvA_bed", required=True, help="RSV-A gene BED file")
    parser.add_argument("--rsvB_bed", required=True, help="RSV-B gene BED file")
    parser.add_argument("--out", required=True, help="Output figure file name (png)")
    args = parser.parse_args()
    
    # Load SARS-CoV-2 to get sample list
    df_sars_wide = load_mosdepth(args.sars)
    sars_samples = df_sars_wide.columns[3:].tolist()
    
    # Compute gene coverage & breadth
    df_sars_cov, df_sars_breadth = compute_gene_coverage(df_sars_wide, load_bed(args.sars_bed), sars_samples)
    df_rsvA_cov, df_rsvA_breadth = compute_gene_coverage(load_mosdepth(args.rsvA), load_bed(args.rsvA_bed), sars_samples)
    df_rsvB_cov, df_rsvB_breadth = compute_gene_coverage(load_mosdepth(args.rsvB), load_bed(args.rsvB_bed), sars_samples)
    
    # Plot multi-panel figure
    fig, axes = plt.subplots(3, 1, figsize=(15,12))
    plot_heatmap_with_breadth(df_sars_cov, df_sars_breadth, "SARS-CoV-2", ax=axes[0])
    plot_heatmap_with_breadth(df_rsvA_cov, df_rsvA_breadth, "RSV-A", ax=axes[1])
    plot_heatmap_with_breadth(df_rsvB_cov, df_rsvB_breadth, "RSV-B", ax=axes[2])
    
    # Add explanatory text at the bottom
    fig.text(0.5, 0.02, "Cells with black borders indicate genes where ≥80% of the gene has ≥5× coverage",
             ha='center', fontsize=10)
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(args.out, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()

