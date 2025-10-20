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
    """Parse read_summary.tsv → list of HA/N serotype strings."""
    if not os.path.exists(file_path):
        return []
    df = pd.read_csv(file_path, sep="\t")
    if df.empty:
        return []

    df["HA"] = df["read_assignment"].str.extract(r"(H\d+)")
    df["NA"] = df["read_assignment"].str.extract(r"(N\d+)")

    serotypes = []
    for _, row in df.iterrows():
        if pd.notna(row["HA"]):
            serotypes.append(row["HA"])
        if pd.notna(row["NA"]):
            serotypes.append(row["NA"])
    return serotypes


def parse_reads_tsv(file_path, influenza=False):
    """Parse summary TSV → dict of {Sample_ID: total reads (sum R1+R2)}."""
    if file_path is None or not os.path.exists(file_path):
        return {}
    df = pd.read_csv(file_path, sep="\t", dtype=str)

    if influenza:
        df = df[df["file"].str.contains("11308")]
        if df.empty:
            raise ValueError(f"No sample rows with '11308' found in {file_path}")

    df["Sample_ID"] = df["file"].str.extract(r"(WW\d+-\d+)")[0]
    df["reads"] = pd.to_numeric(df["num_seqs"], errors="coerce").fillna(0)
    return df.groupby("Sample_ID")["reads"].sum().to_dict()


def generate_heatmap(
    df, all_serotypes, all_dates, site_order, output_prefix, suffix, title
):
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
        hspace=0.4,  # uniform inter-panel spacing
        top=0.95,  # reduce top margin
        bottom=0.05,  # reduce bottom margin
    )

    vmax_global = df["count"].max() if not df.empty else 0

    for idx, site in enumerate(site_order):
        ax = fig.add_subplot(gs[idx, 0])
        site_df = df[df["site"] == site]
        site_pivot = site_df.pivot_table(
            index="date_str", columns="serotype", values="count", fill_value=0
        )
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

        ### 🔧 CHANGED: Improve label readability
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.tick_params(axis="y", labelrotation=0, labelsize=8)

    ### 🔧 CHANGED: Replace tight_layout() with explicit adjustment
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.07, right=0.95)

    out_png = f"{output_prefix}_{suffix}_heatmap.png"
    plt.savefig(out_png, dpi=300)
    plt.close(fig)


def main(
    base_dir,
    metadata_path,
    output_prefix,
    normalize,
    influenza_reads_file,
    total_reads_file,
    filters,
):
    # -------------------------
    # Load and clean metadata
    # -------------------------
    md = Metadata(metadata_path)
    filter_dict = dict(f.split("=") for f in filters) if filters else {}

    # Prefer batch_name from filters, fallback to base_dir
    if "batch_name" in filter_dict:
        batch_name = filter_dict["batch_name"]
    else:
        batch_name = os.path.basename(base_dir.rstrip("/"))

    # Apply metadata filters (excluding batch_name for now)
    meta_df = (
        md.filter(**{k: v for k, v in filter_dict.items() if k != "batch_name"})
        if filters
        else md.df.copy()
    )

    # Restrict to matching batch_name
    meta_df = meta_df[meta_df["batch_name"] == batch_name]
    if meta_df.empty:
        raise ValueError(f"No metadata rows match batch '{batch_name}'")

    # Clean metadata fields
    meta_df["sample_id"] = (
        meta_df["sample_id"].astype(str).str.replace(" ", "-").str.replace("_", "-")
    )
    meta_df["collection_site"] = (
        meta_df["collection_site"]
        .astype(str)
        .str.replace(" ", "-")
        .str.replace("_", "-")
    )
    meta_df["collection_date"] = pd.to_datetime(
        meta_df["collection_date"], errors="coerce"
    )

    # Prepare site/date mapping
    sample_info = (
        meta_df[["sample_id", "collection_site", "collection_date"]]
        .drop_duplicates()
        .rename(columns={"collection_site": "site", "collection_date": "date"})
    ).dropna(subset=["date"])

    # -------------------------
    # Define site order (north → south)
    # -------------------------
    matched_sites = sorted(sample_info["site"].unique())
    prefix_order = ["NHA", "IHA", "MV", "VIHA"]
    ordered_sites = []
    for pref in prefix_order:
        ordered_sites += sorted(
            [s for s in matched_sites if str(s).upper().startswith(pref + "-")]
        )
    ordered_sites += sorted(
        [
            s
            for s in matched_sites
            if not any(str(s).upper().startswith(p + "-") for p in prefix_order)
        ]
    )
    site_order = ordered_sites

    # -------------------------
    # Sort by date then site order
    # -------------------------
    site_rank = {site: i for i, site in enumerate(site_order)}
    sample_info["site_rank"] = (
        sample_info["site"].map(site_rank).fillna(len(site_order))
    )
    sample_info = sample_info.sort_values(["date", "site_rank"]).reset_index(drop=True)

    # -------------------------
    # Collect serotype data
    # -------------------------
    all_serotypes = [f"H{i}" for i in range(1, 17)] + [f"N{i}" for i in range(1, 10)]
    serotype_records = []
    entries = os.listdir(base_dir)

    for _, row in sample_info.iterrows():
        sid, site, date = row["sample_id"], row["site"], row["date"]
        candidates = [f for f in entries if f.startswith(sid)]
        if not candidates:
            candidates = [f for f in entries if sid in f or sid.replace("-", "") in f]
        serotypes = []
        if candidates:
            sample_folder = candidates[0]
            serotype_path = os.path.join(
                base_dir, sample_folder, "11308", "alignment", "serotype_assignment"
            )
            read_summary_file = os.path.join(
                serotype_path, f"{sample_folder}_read_summary.tsv"
            )
            serotypes = parse_read_summary(read_summary_file)
        for ser in serotypes:
            serotype_records.append(
                {"sample": sid, "site": site, "date": date, "serotype": ser, "count": 1}
            )

    # -------------------------
    # Combine baseline + observed
    # -------------------------
    baseline_records = [
        {
            "sample": r.sample_id,
            "site": r.site,
            "date": r.date,
            "serotype": s,
            "count": 0,
        }
        for _, r in sample_info.iterrows()
        for s in all_serotypes
    ]
    df = pd.DataFrame(baseline_records)
    if serotype_records:
        df = pd.concat([df, pd.DataFrame(serotype_records)], ignore_index=True)

    df = df.groupby(["sample", "site", "date", "serotype"], as_index=False)[
        "count"
    ].sum()
    df["date_str"] = df["date"].dt.strftime("%Y-%m-%d")

    # -------------------------
    # Optional normalization
    # -------------------------
    influenza_reads = (
        parse_reads_tsv(influenza_reads_file, influenza=True)
        if normalize and influenza_reads_file
        else {}
    )
    total_reads = (
        parse_reads_tsv(total_reads_file, influenza=False)
        if normalize and total_reads_file
        else {}
    )

    if normalize:
        # Per-influenza CPM
        df_influenza = df.copy()
        df_influenza["denominator"] = df_influenza["sample"].map(
            lambda s: influenza_reads.get(s, 0)
        )
        df_influenza["count"] = df_influenza.apply(
            lambda r: (
                (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0
            ),
            axis=1,
        )
        generate_heatmap(
            df_influenza,
            all_serotypes,
            sorted(df["date"].unique()),
            site_order,
            output_prefix,
            "per_influenza_CPM",
            "Per-influenza CPM",
        )

        # Per-total CPM
        df_total = df.copy()
        df_total["denominator"] = df_total["sample"].map(
            lambda s: total_reads.get(s, 0)
        )
        df_total["count"] = df_total.apply(
            lambda r: (
                (r["count"] / r["denominator"] * 1e6) if r["denominator"] > 0 else 0
            ),
            axis=1,
        )
        generate_heatmap(
            df_total,
            all_serotypes,
            sorted(df["date"].unique()),
            site_order,
            output_prefix,
            "per_total_CPM",
            "Per-total CPM",
        )
    else:
        generate_heatmap(
            df,
            all_serotypes,
            sorted(df["date"].unique()),
            site_order,
            output_prefix,
            "raw_counts",
            "Serotype Count",
        )

    # -------------------------
    # Write summary TSV
    # -------------------------
    summary_rows = []
    for sid, site, date in sample_info[["sample_id", "site", "date"]].values:
        row = {
            "sample_id": sid,
            "collection_date": date.strftime("%Y-%m-%d"),
            "site": site,
        }
        subset = df[(df["sample"] == sid) & (df["site"] == site)]
        for s in all_serotypes:
            row[s] = (
                subset.loc[subset["serotype"] == s, "count"].sum()
                if not subset.empty
                else 0
            )
        row["total_reads"] = total_reads.get(sid, 0)
        row["reads_mapping_influenza"] = influenza_reads.get(sid, 0)
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    out_tsv = f"{output_prefix}_serotype_summary.tsv"
    summary_df.to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Influenza serotype heatmap generator (metadata-integrated, sorted by date and site)"
    )
    parser.add_argument("--base_dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--influenza_reads_summary")
    parser.add_argument("--total_reads_summary")
    parser.add_argument(
        "--filter",
        nargs="+",
        help="Optional metadata filters, e.g. season=2025_2026 extraction_source=research period=in_season batch_name=SAFEGUARD_ww_probe_batch_23_2025_2026",
    )
    args = parser.parse_args()

    main(
        args.base_dir,
        args.metadata,
        args.output_prefix,
        args.normalize,
        args.influenza_reads_summary,
        args.total_reads_summary,
        args.filter,
    )
