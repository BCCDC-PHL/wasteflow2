#!/usr/bin/env python3
import pandas as pd
import re
import argparse
import sys


# ----------------------------
# Sample class
# ----------------------------
class Sample:
    def __init__(self, record):
        self.sample_id = record["sample_id"]
        self.season = record["season"]
        self.period = record["period"]
        self.extraction_source = record["extraction_source"]
        self.batch_name = record["batch_name"]
        self.collection_site = record["collection_site"]
        self.collection_date = record["collection_date"]
        self.sequencing_run = record["sequencing_run"]
        self.test_code = record["test_code"]
        self.results = record["results"]

    def __repr__(self):
        return f"Sample({self.sample_id}, {self.season}, {self.period}, {self.extraction_source})"


# ----------------------------
# Metadata parser
# ----------------------------
class Metadata:
    def __init__(self, metadata_file):
        print(f"Loading metadata file: {metadata_file}")
        self.df = pd.read_csv(metadata_file, sep="\t", dtype=str)

        # --- Clean column names ---
        self.df.columns = (
            self.df.columns.str.strip()
            .str.replace(" ", "_")
            .str.replace("-", "_")
            .str.lower()
        )

        # --- Clean values ---
        self.df = self.df.applymap(lambda x: str(x).strip() if pd.notnull(x) else x)

        # 1. sample_ID → WWXX-XXXX (replace spaces)
        if "sample_id" in self.df.columns:
            self.df["sample_id"] = self.df["sample_id"].str.replace(
                " ", "-", regex=False
            )
        else:
            print(
                "ERROR: 'sample_id' column not found after cleaning.", file=sys.stderr
            )
            sys.exit(1)

        # 2. season should be YYYY or YYYY_YYYY
        self.df["season"] = self.df["season"].apply(
            lambda x: (
                x if pd.isna(x) else (x if re.match(r"^\d{4}(_\d{4})?$", x) else pd.NA)
            )
        )

        # 3. period → off_season / in_season
        self.df["period"] = (
            self.df["period"]
            .str.lower()
            .replace({"offseason": "off_season", "inseason": "in_season"})
        )
        self.df["period"] = self.df["period"].where(
            self.df["period"].isin(["off_season", "in_season"])
        )

        # 4. remove spaces everywhere
        self.df = self.df.applymap(
            lambda x: x.replace(" ", "_") if isinstance(x, str) else x
        )

        # 5. collection_site should be dash-separated
        self.df["collection_site"] = self.df["collection_site"].str.replace("_", "-")

        # 6. convert collection_date from mm/dd/yyyy → yyyy/mm/dd
        self.df["collection_date"] = pd.to_datetime(
            self.df["collection_date"], format="%m/%d/%Y", errors="coerce"
        ).dt.strftime("%Y/%m/%d")

        # --- Convert to Sample objects ---
        self.samples = [Sample(row) for _, row in self.df.iterrows()]

    # ----------------------------
    # Filtering helper
    # ----------------------------
    def filter(self, **kwargs):
        """
        Filter samples by metadata columns, e.g.:
        md.filter(season='2025_2026', period='in_season', extraction_source='research')
        """
        df_filtered = self.df.copy()
        for key, value in kwargs.items():
            if key in df_filtered.columns:
                df_filtered = df_filtered[df_filtered[key] == value]
            else:
                print(f"Warning: '{key}' not found in metadata columns.")
        return df_filtered


# ----------------------------
# Main CLI
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Metadata parser and cleaner for SAFEGUARD workflow."
    )
    parser.add_argument(
        "--metadata_file", required=True, help="Path to metadata TSV file"
    )
    parser.add_argument(
        "--filter",
        nargs="*",
        help="Optional filters, e.g. season=2025_2026 period=in_season extraction_source=research",
    )
    args = parser.parse_args()

    md = Metadata(args.metadata_file)

    # Only print filtered samples if filters are provided
    if args.filter:
        filters = dict(item.split("=") for item in args.filter)
        filtered_df = md.filter(**filters)
        if filtered_df.empty:
            print("⚠️ No samples matched the filter criteria.")
        else:
            print(filtered_df.to_string(index=False))
