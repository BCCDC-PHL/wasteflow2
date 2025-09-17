#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

def parse_read_summary(file_path):
    """Parse read_summary.tsv → list of HA/N serotype strings."""
    if not os.path.exists(file_path):
        return []
    df = pd.read_csv(file_path, sep="\t")
    if df.empty:
        return []

    df["HA"] = df["serotype"].str.extract(r"(H\d+)")
    df["NA"] = df["serotype"].str.extract(r"(N\d+)")

    serotypes = []
    for _, row in df.iterrows():
        if pd.notna(row["HA"]):
            serotypes.append(row["HA"])
        if pd.notna(row["NA"]):
            serotypes.append(row["NA"])
    
    # print(f"🛈 Parsed serotypes from {file_path}: {serotypes}")
    return serotypes

def parse_reads_tsv(file_path, influenza=False):
    """Parse summary TSV → dict of {Sample_ID: total reads (sum R1+R2)}."""
    if file_path is None or not os.path.exists(file_path):
        return {}

    df = pd.read_csv(file_path, sep="\t", dtype=str)

    if influenza:
        # Only keep rows with 11308 (influenza extraction)
        df = df[df["file"].str.contains("11308")]
        if df.empty:
            raise ValueError(f"No sample rows with '11308' found in {file_path}")
    # Extract Sample_ID
    df["Sample_ID"] = df["file"].str.extract(r"(WW\d+-\d+)")[0]
    # Convert num_seqs to numeric
    df["reads"] = pd.to_numeric(df["num_seqs"], errors="coerce").fillna(0)
    # Sum all reads per sample (R1 + R2)
    per_sample_reads = df.groupby("Sample_ID")["reads"].sum().to_dict()
    return per_sample_reads

def generate_heatmap(df, all_serotypes, all_dates, site_order, output_prefix, suffix, title):
    # Create colormap with white as zero
    base_cmap = sns.light_palette("seagreen", as_cmap=True)
    cmap = ListedColormap(['#ffffff'] + [base_cmap(i) for i in range(1, 256)])

    fig_width = max(10, len(all_serotypes) * 0.4)
    fig_height = max(4, len(all_dates) * 0.5 * max(1, len(site_order)))

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(nrows=len(site_order), ncols=1,
                           height_ratios=[len(all_dates)] * len(site_order),
                           hspace=1.0)

    vmax_global = df["count"].max()
    print(f"🛈 Generating heatmap '{title}' with vmax={vmax_global}")

    for idx, site in enumerate(site_order):
        ax = fig.add_subplot(gs[idx, 0])
        site_df = df[df["site"] == site]
        site_pivot = site_df.pivot_table(index="date_str", columns="serotype", values="count", fill_value=0)
        site_pivot = site_pivot.reindex(index=[d.strftime("%Y-%m-%d") for d in all_dates],
                                        columns=all_serotypes, fill_value=0)

        sns.heatmap(site_pivot, ax=ax, cmap=cmap, vmin=0, vmax=vmax_global,
                    linewidths=0.5, linecolor='gray',
                    cbar_kws={"label": title, "format": "%.0f"})

        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_title(site, fontsize=12, loc="left", pad=16)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")

    plt.tight_layout()
    out_png = f"{output_prefix}_{suffix}_heatmap.png"
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"✅ {suffix.capitalize()} heatmap saved -> {out_png} (vmax={int(vmax_global)})")

def main(base_dir, metadata_path, output_prefix, normalize, influenza_reads_file, total_reads_file):
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    metadata = metadata.rename(columns={
        "Patient First Name": "Collection_site",
        "Collection Date": "Collection_date"
    })
    metadata = metadata.dropna(subset=["Collection_site", "Collection_date"])
    metadata["Collection_date"] = pd.to_datetime(metadata["Collection_date"], errors="coerce")
    metadata = metadata.dropna(subset=["Collection_date"])

    batch_name = os.path.basename(base_dir.rstrip("/"))
    metadata = metadata[metadata["Batch_name"] == batch_name]
    if metadata.empty:
        raise ValueError(f"No metadata rows match batch '{batch_name}'")

    metadata["Sample_ID"] = metadata["Sample_ID"].astype(str).str.replace(" ", "-", regex=False)
    sample_info = metadata[["Sample_ID", "Collection_site", "Collection_date"]].drop_duplicates()

    # print(f"\n🛈 Metadata filtered for batch '{batch_name}' ({len(sample_info)} samples):\n", sample_info.head())

    all_serotypes = [f"H{i}" for i in range(1, 17)] + [f"N{i}" for i in range(1, 10)]
    serotype_records = []

    for _, row in sample_info.iterrows():
        sample_id, site, date = row["Sample_ID"], row["Collection_site"], row["Collection_date"]

        candidate_folders = [f for f in os.listdir(base_dir) if f.startswith(sample_id)]
        serotypes = []
        if candidate_folders:
            sample_folder = candidate_folders[0]
            serotype_path = os.path.join(base_dir, sample_folder, "11308", "alignment", "serotype_assignment")
            read_summary_file = os.path.join(serotype_path, f"{sample_folder}_read_summary.tsv")
            serotypes = parse_read_summary(read_summary_file)

        # print(f"🛈 Sample {sample_id} ({site}, {date.date()}) -> serotypes: {serotypes}")

        for ser in serotypes:
            serotype_records.append({"sample": sample_id, "site": site, "date": date, "serotype": ser, "count": 1})

    serotype_df = pd.DataFrame(serotype_records)

    baseline_records = []
    for _, row in sample_info.iterrows():
        for ser in all_serotypes:
            baseline_records.append({
                "sample": row["Sample_ID"],
                "site": row["Collection_site"],
                "date": row["Collection_date"],
                "serotype": ser,
                "count": 0
            })
    baseline_df = pd.DataFrame(baseline_records)

    df = pd.concat([baseline_df, serotype_df], ignore_index=True)
    df = df.groupby(["sample", "site", "date", "serotype"], as_index=False)["count"].sum()
    df["date_str"] = df["date"].dt.strftime("%Y-%m-%d")

    # print(f"\n🛈 Combined DataFrame shape: {df.shape}")
    # print(df.head(10))

    matched_sites = sorted(df["site"].unique())
    all_dates = sorted(df["date"].unique(), reverse=True)
    site_order = (
        sorted([s for s in matched_sites if s.startswith("MV-")]) +
        sorted([s for s in matched_sites if s.startswith("VIHA-")]) +
        sorted([s for s in matched_sites if s.startswith("IHA-")]) +
        sorted([s for s in matched_sites if s.startswith("NHA-")]) +
        sorted([s for s in matched_sites if not (s.startswith("MV-") or s.startswith("VIHA-") or s.startswith("IHA-") or s.startswith("NHA-"))])
    )

    if normalize:
        # Parse reads correctly
        influenza_reads = parse_reads_tsv(influenza_reads_file, influenza=True)
        total_reads = parse_reads_tsv(total_reads_file, influenza=False)

        # 🛈 Debug: show denominators per sample
        print("\n🛈 Influenza reads per sample (denominator for per-influenza CPM):")
        for sample, reads in influenza_reads.items():
            print(f"  {sample}: {reads} reads")

        print("\n🛈 Total reads per sample (denominator for per-total CPM):")
        for sample, reads in total_reads.items():
            print(f"  {sample}: {reads} reads")

        # ---- Per-influenza CPM ----
        df_influenza = df.copy()
        df_influenza["denominator"] = df_influenza["sample"].map(lambda s: influenza_reads.get(s, 0))
        df_influenza["count"] = df_influenza.apply(
            lambda r: (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0, axis=1
        )
        print("\n🛈 Top 10 Per-influenza CPM values (with denominator):")
        print(df_influenza.sort_values("count", ascending=False).head(10)[["sample", "serotype", "denominator", "count"]])
        generate_heatmap(df_influenza, all_serotypes, all_dates, site_order, output_prefix, "per_influenza_CPM", "Per-influenza CPM")

        # ---- Per-total CPM ----
        df_total = df.copy()
        df_total["denominator"] = df_total["sample"].map(lambda s: total_reads.get(s, 0))
        df_total["count"] = df_total.apply(
            lambda r: (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0, axis=1
        )
        print("\n🛈 Top 10 Per-total CPM values (with denominator):")
        print(df_total.sort_values("count", ascending=False).head(10)[["sample", "serotype", "denominator", "count"]])
        generate_heatmap(df_total, all_serotypes, all_dates, site_order, output_prefix, "per_total_CPM", "Per-total CPM")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Influenza serotype heatmap generator")
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--normalize", action="store_true", help="If set, generate normalized heatmaps")
    parser.add_argument("--influenza_reads_summary", help="TSV with per-sample influenza-mapped reads")
    parser.add_argument("--total_reads_summary", help="TSV with per-sample total post-QC reads")
    args = parser.parse_args()
    main(args.base_dir, args.metadata, args.output_prefix, args.normalize, args.influenza_reads_summary, args.total_reads_summary)