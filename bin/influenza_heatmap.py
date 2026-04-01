#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

# Import metadata parser
from metadata_parser import Metadata


def parse_read_summary(file_path):
    """
    Parse read_assignment TSV → return DataFrame with one row per read-pair and
    columns ['pair_id','HA','NA'].

    Handles both 'qname' style (ends with '/1' or '/2') and 'file' style
    (e.g. sample_read_1.fastq.gz / sample_read_2.fastq.gz).
    Deduplicates by pair_id so each pair counts once.
    """
    if not os.path.exists(file_path):
        return pd.DataFrame(columns=["pair_id", "HA", "NA"])

    df = pd.read_csv(file_path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["pair_id", "HA", "NA"])

    # create pair_id depending on available columns
    if "qname" in df.columns:
        df["pair_id"] = df["qname"].str.replace(r"(/1|/2)$", "", regex=True)
    elif "file" in df.columns:
        df["pair_id"] = (
            df["file"]
            .str.replace(r"_read_[12]\.fastq\.gz$", "", regex=True)
            .str.replace(r"_[12]\.fastq\.gz$", "", regex=True)
            .str.replace(r"_R[12]\.fastq\.gz$", "", regex=True)
        )
    else:
        first_col = df.columns[0]
        df["pair_id"] = df[first_col].astype(str).str.replace(r"(/1|/2|_R[12]|_[12])$", "", regex=True)

    df_unique = df.drop_duplicates(subset=["pair_id"]).copy()
    df_unique["HA"] = df_unique.get("read_assignment", pd.Series()).astype(str).str.extract(r"(H\d+)", expand=False)
    df_unique["NA"] = df_unique.get("read_assignment", pd.Series()).astype(str).str.extract(r"(N\d+)", expand=False)

    return df_unique[["pair_id", "HA", "NA"]].reset_index(drop=True)


def parse_target_reads_tsv(file_path):
    """
    Parse combined seqkit summary TSV → return two dicts:
      1. total_target_reads: sum of all R1 read counts per sample
      2. influenza_reads: sum of all R1 read counts per sample that contain '11308'
    
    Works for *_read_1.fastq.gz, *_1.fastp.fastq.gz, etc.
    Sample IDs are extracted using the WW\d+-\d+ pattern.
    """
    if file_path is None or not os.path.exists(file_path):
        return {}, {}

    df = pd.read_csv(file_path, sep="\t", dtype=str)
    if df.empty or "file" not in df.columns or "num_seqs" not in df.columns:
        return {}, {}

    # Keep only R1 reads
    df = df[df["file"].str.contains(r"(_1\.fastp\.fastq\.gz|_read_1\.fastq\.gz)$", regex=True)].copy()
    if df.empty:
        return {}, {}

    # Extract sample IDs (e.g., WW25-1425)
    df["Sample_ID"] = df["file"].str.extract(r"(WW\d+-\d+)")[0]
    df = df.dropna(subset=["Sample_ID"])

    # Convert num_seqs to numeric
    df["num_seqs"] = pd.to_numeric(df["num_seqs"], errors="coerce").fillna(0)

    # total_target_reads → sum all R1 reads per sample
    total_target_reads = df.groupby("Sample_ID")["num_seqs"].sum().to_dict()

    # influenza_reads → subset with '11308' in filename
    df_influenza = df[df["file"].str.contains("11308", na=False)]
    influenza_reads = df_influenza.groupby("Sample_ID")["num_seqs"].sum().to_dict()

    return total_target_reads, influenza_reads


def find_read_summary_files(base_dir):
    """Recursively find all *_read_summary.tsv files in base_dir."""
    read_summary_files = {}
    for root, _, files in os.walk(base_dir):
        for f in files:
            if f.endswith("_read_summary.tsv"):
                file_path = os.path.join(root, f)
                # Extract sample ID (e.g., WW25-1425)
                sample_id = None
                parts = f.split("_")
                for p in parts:
                    if p.startswith("WW") and "-" in p:
                        sample_id = p
                        break
                if sample_id:
                    read_summary_files[sample_id] = file_path
    return read_summary_files


def generate_heatmap(df, all_serotypes, all_dates, site_order, output_prefix, suffix, title):
    """Draw multi-panel heatmap by site."""
    base_cmap = sns.light_palette("seagreen", as_cmap=True)
    cmap = ListedColormap(["#ffffff"] + [base_cmap(i) for i in range(1, 256)])

    panel_height = 2.2 if len(site_order) > 10 else 2.5
    fig_height = len(site_order) * panel_height
    fig_width = max(12, len(all_serotypes) * 0.4)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(
        nrows=len(site_order),
        ncols=1,
        hspace=0.4,
        top=0.95,
        bottom=0.05,
    )

    vmax_global = df["count"].max() if not df.empty else 0

    for idx, site in enumerate(site_order):
        ax = fig.add_subplot(gs[idx, 0])
        site_df = df[df["site"] == site]
        site_pivot = site_df.pivot_table(index="date_str", columns="serotype", values="count", fill_value=0)
        site_pivot = site_pivot.reindex(
            index=[pd.to_datetime(d).strftime("%Y-%m-%d") for d in all_dates],
            columns=all_serotypes,
            fill_value=0,
        )

        sns.heatmap(
            site_pivot,
            ax=ax,
            cmap=cmap,
            vmin=0,
            vmax=vmax_global,
            linewidths=0.5,
            linecolor="gray",
            cbar_kws={"label": title, "format": "%.0f"},
        )

        ax.set_title(site, fontsize=12, loc="left", pad=12)
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.tick_params(axis="y", labelrotation=0, labelsize=8)

    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.07, right=0.95)
    out_png = f"{args.outdir}/{output_prefix}_{suffix}_heatmap.png" if args.outdir else f"{output_prefix}_{suffix}_heatmap.png"
    plt.savefig(out_png, dpi=300)
    plt.close(fig)


def main(base_dir, metadata_path, output_prefix, normalize, target_read_summary, filters, segment_coverage_file=None):
    # -------------------------
    # Load and clean metadata
    # -------------------------
    md = Metadata(metadata_path)
    filter_dict = dict(f.split("=") for f in filters) if filters else {}
    batch_name = filter_dict.get("batch_name", os.path.basename(base_dir.rstrip("/")))

    meta_df = (
        md.filter(**{k: v for k, v in filter_dict.items() if k != "batch_name"})
        if filters
        else md.df.copy()
    )
    meta_df = meta_df[meta_df["batch_name"] == batch_name]
    if meta_df.empty:
        raise ValueError(f"No metadata rows match batch '{batch_name}'")

    meta_df["sample_id"] = meta_df["sample_id"].astype(str).str.replace(" ", "-").str.replace("_", "-")
    meta_df["collection_site"] = meta_df["collection_site"].astype(str).str.replace(" ", "-").str.replace("_", "-")
    meta_df["collection_date"] = pd.to_datetime(meta_df["collection_date"], errors="coerce")

    sample_info = (
        meta_df[["sample_id", "collection_site", "collection_date"]]
        .drop_duplicates()
        .rename(columns={"collection_site": "site", "collection_date": "date"})
    ).dropna(subset=["date"])

    # -------------------------
    # Define site order
    # -------------------------
    matched_sites = sorted(sample_info["site"].unique())
    prefix_order = ["NHA", "IHA", "MV", "VIHA"]
    ordered_sites = []
    for pref in prefix_order:
        ordered_sites += sorted([s for s in matched_sites if str(s).upper().startswith(pref + "-")])
    ordered_sites += sorted([s for s in matched_sites if not any(str(s).upper().startswith(p + "-") for p in prefix_order)])
    site_order = ordered_sites

    site_rank = {site: i for i, site in enumerate(site_order)}
    sample_info["site_rank"] = sample_info["site"].map(site_rank).fillna(len(site_order))
    sample_info = sample_info.sort_values(["date", "site_rank"]).reset_index(drop=True)

    # -------------------------
    # Find summaries and serotype files
    # -------------------------
    summary_file_map = find_read_summary_files(base_dir)

    all_serotypes = [f"H{i}" for i in range(1, 17)] + [f"N{i}" for i in range(1, 10)]
    serotype_records = []

    for _, row in sample_info.iterrows():
        sid, site, date = row["sample_id"], row["site"], row["date"]
        df_pairs = pd.DataFrame(columns=["pair_id", "HA", "NA"])

        for key, file_path in summary_file_map.items():
            if sid in key or key in sid or sid.replace("-", "") in key:
                df_pairs = parse_read_summary(file_path)
                break

        if not df_pairs.empty:
            ha_counts = df_pairs["HA"].value_counts(dropna=True).to_dict()
            na_counts = df_pairs["NA"].value_counts(dropna=True).to_dict()

            for ser, cnt in ha_counts.items():
                if pd.notna(ser):
                    serotype_records.append({"sample": sid, "site": site, "date": date, "serotype": ser, "count": int(cnt)})

            for ser, cnt in na_counts.items():
                if pd.notna(ser):
                    serotype_records.append({"sample": sid, "site": site, "date": date, "serotype": ser, "count": int(cnt)})

    if segment_coverage_file and os.path.exists(segment_coverage_file):
        segment_df = pd.read_csv(segment_coverage_file, sep="\t")
        ha_df = segment_df[(segment_df["Segment"] == "HA") & (segment_df["percent_covered"] > 0)]
        valid_ha_samples = set(zip(ha_df["Sample"], ha_df["Pathogen"]))
        ha_to_pathogen = {"H1": "H1N1", "H3": "H3N2", "H5": "H5N1"}

        filtered_records = []
        for rec in serotype_records:
            if rec["serotype"] in ha_to_pathogen:
                pathogen = ha_to_pathogen[rec["serotype"]]
                if (rec["sample"], pathogen) in valid_ha_samples:
                    filtered_records.append(rec)
            else:
                filtered_records.append(rec)
        serotype_records = filtered_records
    
    # -------------------------
    # Combine + normalize + plot
    # -------------------------
    baseline_records = [
        {"sample": r.sample_id, "site": r.site, "date": r.date, "serotype": s, "count": 0}
        for _, r in sample_info.iterrows()
        for s in all_serotypes
    ]
    df = pd.DataFrame(baseline_records)
    if serotype_records:
        df = pd.concat([df, pd.DataFrame(serotype_records)], ignore_index=True)

    df = df.groupby(["sample", "site", "date", "serotype"], as_index=False)["count"].sum()
    df["date_str"] = df["date"].dt.strftime("%Y-%m-%d")

    # -------------------------
    # Parse target read summary
    # -------------------------
    total_target_reads, influenza_reads = parse_target_reads_tsv(target_read_summary) if normalize and target_read_summary else ({}, {})

    if normalize:
        df_influenza = df.copy()
        df_influenza["denominator"] = df_influenza["sample"].map(lambda s: influenza_reads.get(s, 0))
        df_influenza["count"] = df_influenza.apply(lambda r: (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0, axis=1)
        generate_heatmap(df_influenza, all_serotypes, sorted(df["date"].unique()), site_order, output_prefix, "per_influenza_CPM", "Per-influenza CPM")

        df_total = df.copy()
        df_total["denominator"] = df_total["sample"].map(lambda s: total_target_reads.get(s, 0))
        df_total["count"] = df_total.apply(lambda r: (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0, axis=1)
        generate_heatmap(df_total, all_serotypes, sorted(df["date"].unique()), site_order, output_prefix, "per_total_CPM", "Per-total CPM")
    else:
        generate_heatmap(df, all_serotypes, sorted(df["date"].unique()), site_order, output_prefix, "raw_counts", "Serotype Count")

    # -------------------------
    # Write summary TSV
    # -------------------------
    summary_rows = []
    for sid, site, date in sample_info[["sample_id", "site", "date"]].values:
        row = {"sample_id": sid, "collection_date": date.strftime("%Y-%m-%d"), "site": site}
        subset = df[(df["sample"] == sid) & (df["site"] == site)]
        for s in all_serotypes:
            row[s] = int(subset.loc[subset["serotype"] == s, "count"].sum()) if not subset.empty else 0
        row["total_target_reads"] = total_target_reads.get(sid, 0)
        row["reads_mapping_influenza"] = influenza_reads.get(sid, 0)
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    out_tsv = f"{args.outdir}/{output_prefix}_serotype_summary.tsv" if args.outdir else f"{output_prefix}_serotype_summary.tsv"
    summary_df.to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Influenza serotype heatmap generator (metadata-integrated)")
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=False, help="Optional output directory (defaults to current dir)")
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--target_read_summary", required=False, help="Single seqkit summary TSV containing reads for all targets")
    parser.add_argument("--segment_coverage_file", required=False, help="TSV with percent coverage per segment for each pathogen") 
    parser.add_argument("--filter", nargs="+", help="Optional metadata filters, e.g. season=2025_2026 extraction_source=research period=in_season batch_name=SAFEGUARD_ww_probe_batch_23_2025_2026")
    args = parser.parse_args()

    main(args.base_dir, args.metadata, args.output_prefix, args.normalize, args.target_read_summary, args.filter, args.segment_coverage_file)