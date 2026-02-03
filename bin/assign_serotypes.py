#!/usr/bin/env python3

import os
import re
import sys
import logging
import argparse
from typing import Dict, List, Tuple

import pysam
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger("iavs_logger")
logging.basicConfig(level=logging.INFO)


def bam_to_score(in_bam: str, chunk_size: int = 100_000) -> pl.DataFrame:
    dfs = []
    buffer = []

    schema = [
        'qname', 'rname',
        'read_length', 'aln_prop',
        'aln_len', 'aln_acc',
        'orientation', 'pair_name'
    ]

    with pysam.AlignmentFile(in_bam, 'r') as ministream:
        for record in ministream:
            readLength = record.infer_read_length()
            alignLength = record.query_alignment_length

            if record.is_paired:
                orient = "R1" if record.is_read1 else "R2"
            else:
                orient = "N"

            endtup = (".1", "/1", ".2", "/2")
            if record.query_name.endswith(endtup):
                pairn = record.query_name[:-2]
            else:
                pairn = record.query_name.strip()

            try:
                alignEdit = dict(record.get_tags()).get('NM', 0)
            except Exception:
                alignEdit = 0

            alignProp = alignLength / readLength if readLength else 0
            alignAcc = (alignLength - alignEdit) / alignLength if alignLength else 0

            if alignLength >= 100 and alignProp >= 0.9 and alignAcc >= 0.8:
                buffer.append([
                    record.query_name,
                    record.reference_name,
                    readLength,
                    alignProp,
                    alignLength,
                    alignAcc,
                    orient,
                    pairn
                ])

            # 🚀 Flush chunk
            if len(buffer) >= chunk_size:
                dfs.append(
                    pl.DataFrame(buffer, schema=schema, orient="row")
                )
                buffer.clear()

    # flush remainder
    if buffer:
        dfs.append(
            pl.DataFrame(buffer, schema=schema, orient="row")
        )

    if not dfs:
        return pl.DataFrame(schema=schema)

    return pl.concat(dfs)


def load_flu_info(db_info_path: str) -> pl.DataFrame:
    df = pl.read_csv(
        db_info_path, 
        separator='\t', 
        has_header=True,
        schema_overrides={"accession": pl.Utf8, "serotype": pl.Utf8, "segment": pl.Utf8}
    )
    # Ensure required columns
    needed = {"accession", "serotype", "segment"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Flu info file missing columns: {missing}")
    return df.select(["accession", "serotype", "segment"]).with_columns(
        pl.col("accession").cast(pl.Utf8),
        pl.col("serotype").cast(pl.Utf8),
        pl.col("segment").cast(pl.Utf8),
    )


def pair_alns(alndf: pl.DataFrame) -> pl.DataFrame:
    mini2_p_df = alndf.with_columns(
        (pl.col("aln_acc") * pl.col("aln_prop")).alias("aln_score")
    ).sort(
        ['qname', 'orientation', 'aln_score'], 
        descending = True
    ).group_by(
        ['qname', 'pair_name', 'orientation', 'rname']
    ).agg(
        pl.col('read_length').first().alias('read_length'),
        pl.col('aln_prop').first().alias('aln_prop'),
        pl.col('aln_len').first().alias('aln_len'),
        pl.col('aln_acc').first().alias('aln_acc'),
        pl.col('aln_score').first().alias('aln_score')
    ).group_by(
        ['pair_name', 'rname']
    ).agg(
        pl.col('read_length').sum().alias('pair_length'),
        pl.col('read_length').count().alias('pair_count'),
        (pl.col("aln_acc").mean() * pl.col('aln_prop').mean() * pl.col('read_length').sum()).alias('pair_score'),  
    )
    return mini2_p_df


def compute_assignment(df_align: pl.DataFrame, flu_info: pl.DataFrame, score_thresh: float = 90) -> pl.DataFrame:
    if df_align.is_empty():
        return df_align

    merged = df_align.join(flu_info, left_on="rname", right_on="accession", how="inner")

    try:
        listed_clf_df = merged.sort(
            ['pair_name', 'pair_score', 'pair_length', 'serotype'], 
            descending = True
        ).group_by(['pair_name', 'serotype']).agg(
            pl.col('pair_score').first().alias('sero_score'),
            pl.col('pair_length').first().alias('sero_al_length'),
            pl.col('pair_count').first().alias('sero_read_count')
        ).sort(
            ['pair_name', 'sero_score'], 
            descending = True
        ).group_by(['pair_name']).agg(
            pl.col('sero_score').first().alias('first_score'),
            pl.col('sero_score').slice(1,1).first().alias('second_score'),
            pl.col('serotype').first().alias('first_sero'),
            pl.col('serotype').slice(1,1).first().alias('second_sero'),
            pl.col('sero_read_count').first().alias('first_nreads'),
            pl.col('sero_read_count').slice(1,1).first().alias('second_nreads'),
            pl.col('sero_al_length').first().alias('first_alength'),
            pl.col('sero_al_length').slice(1,1).first().alias('second_alength'),
        ).with_columns([
            # decision logic: best serotype has to be better than any other serotype
            pl.when(
                (
                    (pl.col("first_score") - 1 >= pl.col("second_score")) | \
                    pl.col("second_score").is_null() | \
                    (pl.col('first_sero') == pl.col('second_sero'))
                )
            ).then(
                pl.col("first_sero")
            ).otherwise(
                pl.lit("ambiguous")
            ).alias("read_assignment")
        ]).filter(
            pl.col("first_score") >= float(score_thresh)
        )
        return listed_clf_df

    except Exception as e:
        logger.warning(f"could not parse taxonomy from alignments")
        logger.warning(e)
        return pl.DataFrame()


def write_outputs(sum_df: pl.DataFrame, sample: str, out_dir: str) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    if sum_df.is_empty():
        return outputs

    summary_path = os.path.join(out_dir, f"{sample}_per_read_summary.tsv")
    sum_df.write_csv(summary_path, separator='\t', include_header=True)
    outputs["summary_tsv"] = summary_path

    # Bar plot of read_assignment (no pandas/pyarrow)
    counts_df = sum_df.group_by("read_assignment").agg(
        pl.col("first_nreads").sum().alias("read_count")
    ).sort(['read_assignment'])
    labels = counts_df["read_assignment"].to_list()
    values = counts_df["read_count"].to_list()
    plt.figure(figsize=(8, 4))
    plt.bar(labels, values, color="#4C78A8")
    plt.xticks(rotation=90)
    plt.ylabel("read count")
    plt.tight_layout()
    plot_path = os.path.join(out_dir, f"{sample}_read_serotype_assignment.pdf")
    plt.savefig(plot_path)
    plt.close()
    outputs["assignment_plot_pdf"] = plot_path

    sero_path = os.path.join(out_dir, f"{sample}_per_serotype_summary.tsv")
    counts_df.write_csv(sero_path, separator='\t', include_header=True)
    outputs["sero_tsv"] = sero_path

    # Per-serotype read lists: one row per qname already
    best_per_read = sum_df

    # Write one file per assignment label
    parts = best_per_read.partition_by("read_assignment", as_dict=True, maintain_order=True)
    for label, g in parts.items():
        label = str(label[0])
        rname_path = os.path.join(out_dir, f"{sample}_{label}.txt")
        pl.DataFrame({
            "pair_name": g["pair_name"],
        }).write_csv(rname_path, separator='\t', include_header=False)
        outputs[f"reads_{label}"] = rname_path

    return outputs


def assign_serotypes(bam_path: str, db_info_path: str, sample: str, out_dir: str, score_thresh: float = 90) -> Dict[str, str]:
    """Main function to assign serotypes from BAM alignments"""
    score_df = bam_to_score(bam_path)
    df_align = pair_alns(score_df)
    
    if df_align.is_empty():
        logger.info("No alignments found; skipping assignment outputs.")
        return {}

    flu_info = load_flu_info(db_info_path)
    sum_df = compute_assignment(df_align, flu_info, score_thresh)

    outputs = write_outputs(sum_df, sample, out_dir)
    return outputs


def main():
    parser = argparse.ArgumentParser(
        description='Assign influenza serotypes from BAM alignments'
    )
    parser.add_argument(
        '--bam',
        required=True,
        help='Input BAM file with alignments'
    )
    parser.add_argument(
        '--db-info',
        required=True,
        help='Database info TSV file (accession, serotype, segment columns)'
    )
    parser.add_argument(
        '--sample',
        required=True,
        help='Sample name for output files'
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output directory (default: current directory)'
    )
    parser.add_argument(
        '--score-thresh',
        type=float,
        default=90.0,
        help='Minimum score threshold for assignment (default: 90.0)'
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Run assignment
    logger.info(f"Processing BAM file: {args.bam}")
    logger.info(f"Using database info: {args.db_info}")
    logger.info(f"Sample name: {args.sample}")
    logger.info(f"Score threshold: {args.score_thresh}")
    
    outputs = assign_serotypes(
        bam_path=args.bam,
        db_info_path=args.db_info,
        sample=args.sample,
        out_dir=args.outdir,
        score_thresh=args.score_thresh
    )
    
    if outputs:
        logger.info(f"Generated {len(outputs)} output files:")
        for key, path in outputs.items():
            logger.info(f"  {key}: {path}")
    else:
        logger.warning("No outputs generated - check alignment quality and thresholds")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())