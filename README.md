# BCCDC-PHL/wasteflow2

> **WasteFlow 2.0** — A scalable Nextflow-based workflow for probe-enriched multi-pathogen surveillance in complex environmental samples.

## Introduction
 
**WasteFlow 2.0** is a bioinformatics workflow designed for simultaneous genomic surveillance of multiple respiratory pathogens from wastewater and other complex environmental matrices. Built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html), the pipeline processes probe-enriched metagenomic sequencing data and routes reads through shared core processes and pathogen-specific analytical subworkflows in parallel.
 
The workflow currently supports six viral targets:
 
| Target | Subtypes / Lineages |
|---|---|
| SARS-CoV-2 | All lineages (Pango nomenclature) |
| Influenza A | All subtypes (e.g. H1N1, H3N2, **H5Nx**) |
| Influenza B | Victoria lineages |
| RSV A | Genotype A - A.D |
| RSV B | Genotype B - B.D.E |
| Measles | All genotypes (A–H) |

The modular Nextflow branching architecture means additional viral targets can be incorporated with minimal changes to the core pipeline — simply define a new target configuration and branch entry point (see [Adding New Viral Targets](#adding-new-viral-targets)).
 
The pipeline is built on [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html), and uses Docker, Singularity, or Conda containers to ensure reproducibility. Where possible, modules follow [nf-core/modules](https://github.com/nf-core/modules) conventions.

---
 
## Pipeline Summary
 
WasteFlow 2.0 processes all targets through a common set of core steps before branching into pathogen-specific subworkflows. This design ensures consistent QC and preprocessing across targets while allowing each virus's unique analytical requirements to be handled independently and simultaneously.
 
### Workflow Overview
 
```
                     ┌──────────────────────────────┐
                     │      INPUT: Raw FASTQ Reads   │
                     └─────────────┬────────────────┘
                                   │
                     ┌─────────────▼────────────────┐
                     │   CORE: Read QC & Trimming    │
                     │       (FastQC, fastp)         │
                     └─────────────┬────────────────┘
                                   │
                     ┌─────────────▼────────────────┐
                     │  CORE: Taxonomic Classification│
                     │   (Kraken2 → Bracken → Krona) │
                     └─────────────┬────────────────┘
                                   │
                     ┌─────────────▼────────────────┐
                     │   CORE: Target Read Extraction │
                     │     (per-virus read binning)  │
                     └──┬──────┬──────┬──────┬───┬──┘
                        │      │      │      │   │
          ┌─────────────┘  ┌───┘  ┌───┘  ┌──┘   └──────────┐
          │                │      │      │                   │
   ┌──────▼──────┐  ┌──────▼─┐  ┌─▼─────▼──┐  ┌────────────▼──┐
   │  SARS-CoV-2 │  │ Flu A  │  │  RSV A+B  │  │    Measles    │
   │             │  │        │  │ (→ split  │  │               │
   │             │  │Serotype│  │  A and B) │  │               │
   │             │  │ calling│  │           │  │               │
   └──────┬──────┘  └───┬────┘  └─────┬────┘  └───────┬───────┘
          │             │             │               │
          │         ┌───▼────┐        │               │      ┌────────┐
          │         │H1N1/   │        │               │      │ Flu B  │
          │         │H3N2/   │        │               │      │(Victoria
          │         │H5Nx    │        │               │      │/Yamag.)│
          │         └───┬────┘        │               │      └───┬────┘
          │             │             │               │          │
          └─────────────┴─────────────┴───────────────┴──────────┘
                                      │
                     ┌────────────────▼─────────────────┐
                     │   CORE: Alignment (Minimap2)      │
                     │   BAM processing (SAMtools)       │
                     │   [Flu A/B: segment splitting]    │
                     └──────────────┬───────────────────┘
                                    │
               ┌────────────────────┼──────────────────────┐
               │                    │                       │
  ┌────────────▼──────────┐  ┌──────▼──────────┐  ┌────────▼────────────┐
  │  Coverage Analysis    │  │ Variant Calling  │  │ Lineage Deconvolution│
  │    (Mosdepth)         │  │    (iVar)        │  │  [Flu A/B: HA-only] │
  │                       │  │      ↓           │  │     (Freyja)        │
  │                       │  │ Variant Annotation│  │                     │
  │                       │  │    (snpEff)      │  │                     │
  └────────────┬──────────┘  └──────┬──────────┘  └────────┬────────────┘
               │                    │                       │
               └────────────────────┴───────────────────────┘
                                    │
                     ┌──────────────▼───────────────────┐
                     │   Downstream Analysis             │
                     │   (VIRUS-MVP, custom integrations)│
                     └──────────────┬───────────────────┘
                                    │
                     ┌──────────────▼───────────────────┐
                     │   REPORTING: MultiQC +            │
                     │   Surveillance Summary            │
                     └──────────────────────────────────┘
```
 
### Core Processes (All Targets)
 
These steps run for every sample regardless of viral target:
 
1. **Read QC** — Quality assessment with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. **Adapter & Quality Trimming** — [fastp](https://github.com/OpenGene/fastp)
3. **Taxonomic Classification** — Community-level profiling with [Kraken2](http://ccb.jhu.edu/software/kraken2/), abundance re-estimation with [Bracken](https://ccb.jhu.edu/software/bracken/), and interactive visualisation with [Krona](https://github.com/marbl/Krona)
4. **Target Read Extraction** — Per-virus read binning from classified output; Influenza A reads undergo additional serotype mapping (H1N1 / H3N2 / H5Nx) and Influenza B reads are assigned to Victoria or Yamagata lineage before alignment; RSV A+B reads are extracted jointly then split into RSV-A and RSV-B
5. **Reference Alignment** — Read alignment to per-target reference genomes with [Minimap2](https://github.com/lh3/minimap2)
6. **BAM Processing** — Sorting, indexing, and alignment statistics with [SAMtools](http://www.htslib.org/) (flagstat, idxstats, stats); Influenza A and B BAMs are additionally split by genome segment prior to downstream analysis
7. **Coverage Analysis** — Per-base and per-region depth with [Mosdepth](https://github.com/brentp/mosdepth)
8. **Variant Calling** — Low-frequency variant detection with [iVar](https://github.com/andersen-lab/ivar)
9. **Variant Annotation** — Functional annotation of called variants with [snpEff](https://pcingola.github.io/SnpEff/)
10. **VCF-to-GVF Conversion** — Custom WasteFlow module that converts annotated VCF files to [GVF (Genome Variation Format)](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md) and incorporates mature peptide annotation for biologically precise interpretation of variant impacts
11. **Lineage Deconvolution** — Relative abundance estimation of co-circulating lineages/genotypes with [Freyja](https://github.com/andersen-lab/Freyja) (variants → demix → bootstrap → aggregate); for Influenza A and B, Freyja is run on the HA segment only
12. **Aggregate QC Report** — [MultiQC](https://multiqc.info/) summary across all targets and samples

### Pathogen-Specific Subworkflows
 
The branching architecture means each target's reads are processed independently and simultaneously after extraction. The table below summarises the virus-specific steps that sit alongside the shared alignment, variant calling, and Freyja deconvolution steps:
 
| Target | Pre-alignment steps | Lineage / Genotype tool | Notes |
|---|---|---|---|
| **SARS-CoV-2** | Direct extraction → alignment | [Freyja](https://github.com/andersen-lab/Freyja) | Full-genome Freyja demixing; Pango lineage abundance |
| **Influenza A** | Serotype mapping (H/N database) → H1N1 / H3N2 / H5Nx extraction → alignment → segment splitting | [Freyja](https://github.com/andersen-lab/Freyja) (HA segment) | HA segment extracted post-alignment for lineage deconvolution |
| **Influenza B** | Subtype mapping (Victoria / Yamagata) → extraction → alignment → segment splitting | [Freyja](https://github.com/andersen-lab/Freyja) (HA segment) | HA segment extracted post-alignment for lineage deconvolution |
| **RSV A** | Joint RSV A+B extraction → split RSV-A → alignment | [Freyja](https://github.com/andersen-lab/Freyja) | Full-genome Freyja demixing |
| **RSV B** | Joint RSV A+B extraction → split RSV-B → alignment | [Freyja](https://github.com/andersen-lab/Freyja) | Full-genome Freyja demixing |
| **Measles** | Direct extraction → alignment | [Freyja](https://github.com/andersen-lab/Freyja) | Full-genome Freyja demixing |
 
All targets share the same post-alignment mutation analysis stack: **iVar** (variant calling) → **snpEff** (functional annotation) → **VCF-to-GVF conversion with mature peptide annotation** (custom module) → downstream tools (see [Downstream Analysis](#downstream-analysis)).
 
---
 
## Nextflow Branching Architecture
 
WasteFlow 2.0 leverages Nextflow DSL2's `branch` operator to route reads into per-target subworkflows **simultaneously**, enabling all six pathogens to be processed in a single workflow run without redundant preprocessing. Because Nextflow executes these branches asynchronously, reads assigned to SARS-CoV-2 can be in Freyja lineage demixing at the same time Influenza A reads are being serotyped and RSV reads are being split — all from the same sequencing run.
 
After taxonomic classification and target read extraction, a channel of per-sample, per-target reads is branched as follows:
 
---
 
## Downstream Analysis
 
A core design principle of WasteFlow 2.0 is **interoperability**. The pipeline is deliberately built to produce standardised output formats that integrate directly with tools across the viral genomics ecosystem — without requiring reformatting, custom glue scripts, or manual intervention between steps.
 
### Mutation Analysis and GVF Output
 
The full post-alignment mutation analysis chain for all targets is:
 
```
iVar (variant calling)
  → snpEff (functional annotation)
    → VCF-to-GVF conversion (custom WasteFlow module)
      → VIRUS-MVP / downstream tools
```
 
After variant calling with iVar and functional annotation with snpEff, WasteFlow 2.0 converts annotated VCF files to **GVF (Genome Variation Format)** using a custom built-in module. Critically, this conversion step also incorporates **mature peptide annotation** — an important layer of biological context that standard VCF-based tools do not capture. For viruses such as SARS-CoV-2, Influenza, and RSV, the viral polyprotein is cleaved into distinct mature peptides with specific functional roles; annotating variants at this level enables far more precise interpretation of their potential effects on infectivity, immune evasion, and drug or vaccine target sites.
 
GVF is a well-established community standard in viral genomics, and producing GVF output ensures WasteFlow results are compatible with a broad range of downstream analysis and visualisation platforms.
 
### VIRUS-MVP Integration
 
WasteFlow 2.0's GVF output is designed for direct use with [**VIRUS-MVP**](https://github.com/cidgoh/VIRUS-MVP), an interactive, portable platform for the comprehensive surveillance of viral mutations.
 
VIRUS-MVP links viral mutations to functional annotations, providing insights into their predicted effects on viral infectivity, immune evasion, and protein functionality. It features an interactive interface for visualising mutation distributions, a modular and reproducible genomics workflow, and a curated annotation resource capturing known impacts on viral proteins and host interactions. Users can also import custom functional annotations to tailor analyses to specific research needs or emerging pathogens. Currently supporting SARS-CoV-2, mpox, and expanding to Influenza and RSV — targets that align directly with WasteFlow's surveillance scope.
 
Together, the two platforms form a complete wastewater genomic surveillance and interpretation stack:
 
- **WasteFlow 2.0** answers *what mutations are present* in the community, at what frequency, and in which co-circulating lineages
- **VIRUS-MVP** answers *what those mutations mean* — their functional consequences, public health relevance, and biological context
Both platforms are open-source and developed collaboratively with public health and academic partners.
 
---
 
## Adding New Viral Targets
 
WasteFlow 2.0 is designed to be extensible. Adding a new viral target requires:
 
1. **Reference genome** — Add a reference FASTA and annotation to `assets/references/`
2. **Probe panel BED file** — Add the target capture probe coordinates to `assets/probes/`
3. **Kraken2 / classification database** — Include the target in the classification database or add a custom k-mer set
4. **Branch entry** — Add a new branch condition in `workflows/wasteflow.nf`
5. **Subworkflow** — Create a new subworkflow in `subworkflows/local/<target>/` following the existing template
6. **Configuration** — Add target-specific parameters to `conf/targets.config`
A template for new target subworkflows is provided at `subworkflows/local/template_target/`.
 
---
 
## Quick Start
 
### 1. Install Nextflow
 
```bash
# Requires Java 17 or later
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```
 
### 2. Install a container engine
 
Choose one of: [Docker](https://docs.docker.com/get-docker/), [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps), or [Conda](https://docs.conda.io/en/latest/miniconda.html).
 
### 3. Prepare your samplesheet
 
Create a CSV samplesheet (`samplesheet.csv`) with the following columns:
 
```csv
sample,fastq_1,fastq_2,probe_panel
SAMPLE_01,/path/to/sample01_R1.fastq.gz,/path/to/sample01_R2.fastq.gz,respiratory_v2
SAMPLE_02,/path/to/sample02_R1.fastq.gz,/path/to/sample02_R2.fastq.gz,respiratory_v2
```
 
The `probe_panel` field maps to a panel definition in `conf/panels.config`. Use `respiratory_v2` for the default panel covering all six targets.
 
### 4. Run the pipeline
 
```bash
nextflow run your-org/wasteflow \
    -profile singularity \
    --input samplesheet.csv \
    --outdir ./results 
```
 
To target specific pathogens only (e.g. skip Measles and RSV B):
 
```bash
nextflow run your-org/wasteflow \
    -profile docker \
    --input samplesheet.csv \
    --outdir ./results \
    --targets "sarscov2,influenza_a,influenza_b,rsv_a"
```
 
---
 
## Parameters
 
### Required
 
| Parameter | Description |
|---|---|
| `--input` | Path to samplesheet CSV (see format above) |
| `--outdir` | Output directory for results |
 
### Core Options
 
| Parameter | Default | Description |
|---|---|---|
| `--targets` | `all` | Comma-separated list of targets to run. Options: `sarscov2`, `influenza_a`, `influenza_b`, `rsv_a`, `rsv_b`, `measles`, or `all` |
| `--genome` | `GRCh38` | Host reference genome for read removal |
| `--probe_panel` | `respiratory_v2` | Probe panel to use for target classification |
| `--platform` | `illumina` | Sequencing platform (`illumina` or `nanopore`) |
| `--min_coverage` | `10` | Minimum read depth threshold for consensus calling |
| `--min_allele_freq` | `0.03` | Minimum allele frequency for variant calling |
 
### Skip Options
 
| Parameter | Description |
|---|---|
| `--skip_fastqc` | Skip FastQC read QC |
| `--skip_kraken2` | Skip Kraken2 classification (use if pre-classified) |
| `--skip_snpeff` | Skip snpEff variant annotation |
| `--skip_freyja` | Skip Freyja lineage deconvolution |
| `--skip_multiqc` | Skip MultiQC report generation |
 
### Target-Specific Options
 
| Parameter | Default | Description |
|---|---|---|
| `--freyja_barcodes` | `latest` | Freyja barcode file version or path to custom barcodes |
| `--flu_a_serotypes` | `h1n1,h3n2,h5nx` | Comma-separated Influenza A serotypes to extract and process |
| `--flu_b_lineages` | `victoria,yamagata` | Influenza B lineages to extract and process |
| `--min_coverage` | `10` | Minimum read depth threshold for consensus and variant calling |
| `--min_allele_freq` | `0.03` | Minimum allele frequency for iVar variant calling |
 
Full parameter documentation is available in the [parameter docs](docs/parameters.md).
 
---
 
## Output
 
Results are written to `--outdir` with the following structure:
 
```
results/
├── sarscov2/
│   ├── alignment/          # BAM files and SAMtools stats
│   ├── coverage/           # Mosdepth per-base and summary files
│   ├── variants/           # iVar VCF files
│   ├── annotation/         # snpEff annotated VCFs
│   ├── gvf/                # GVF files with mature peptide annotation
│   └── freyja/             # Demix, bootstrap, and aggregate outputs
├── influenza_a/
│   ├── alignment/          # BAM files per serotype (H1N1, H3N2, H5Nx)
│   ├── segments/           # Per-segment BAMs
│   ├── coverage/           # Mosdepth per-segment coverage
│   ├── variants/           # iVar VCFs per segment
│   ├── annotation/         # snpEff annotated VCFs
│   ├── gvf/                # GVF files with mature peptide annotation
│   └── freyja/             # HA-segment lineage deconvolution
├── influenza_b/
│   ├── alignment/          # BAM files per lineage (Victoria, Yamagata)
│   ├── segments/
│   ├── coverage/
│   ├── variants/
│   ├── annotation/
│   ├── gvf/
│   └── freyja/             # HA-segment lineage deconvolution
├── rsv_a/
│   ├── alignment/
│   ├── coverage/
│   ├── variants/
│   ├── annotation/
│   ├── gvf/
│   └── freyja/
├── rsv_b/
│   ├── alignment/
│   ├── coverage/
│   ├── variants/
│   ├── annotation/
│   ├── gvf/
│   └── freyja/
├── measles/
│   ├── alignment/
│   ├── coverage/
│   ├── variants/
│   ├── annotation/
│   ├── gvf/
│   └── freyja/
├── taxonomy/
│   ├── kraken2/            # Per-sample Kraken2 reports
│   ├── bracken/            # Bracken abundance estimates
│   └── krona/              # Interactive Krona HTML charts
├── multiqc/
│   └── multiqc_report.html # Aggregate QC report
└── pipeline_info/          # Nextflow execution reports and logs
```
 
---
 
## Supported Platforms
 
| Platform | Status | Notes |
|---|---|---|
| Illumina (short-read) | ✅ Supported | Paired-end, 150 bp recommended |
| Nanopore (long-read) | ✅ Supported | Min-pass filtered reads |
| Ion Torrent | 🔶 Experimental | Single-end mode |
 
---
 
## Requirements
 
- **Nextflow** ≥ 24.04.0
- **Java** ≥ 17
- **Container engine**: Docker, Singularity/Apptainer, or Conda
- **Memory**: ≥ 16 GB RAM recommended for multi-target runs
- **Storage**: ≥ 50 GB scratch space per batch of 96 samples
---
 
## Installation
 
```bash
# Clone the repository
git clone https://github.com/your-org/wasteflow.git
cd wasteflow
 
# Test the installation with a minimal dataset
nextflow run . -profile test,docker --outdir ./test_results
```
 
A minimal test dataset covering all six targets is provided in `assets/test_data/`.
 
---
 
## Resource Configuration
 
WasteFlow 2.0 ships with sensible defaults for local workstations, HPC clusters (SLURM, PBS, LSF), and cloud environments (AWS, Google Cloud, Azure). Override with `-profile`:
 
```bash
# Local with Docker
-profile docker
 
# HPC with Singularity + SLURM
-profile singularity,slurm
 
# AWS Batch
-profile aws
 
# Custom — edit conf/custom.config and pass:
-c conf/custom.config
```
 
Institution-specific config profiles can be added to the `conf/` directory or contributed to [nf-core/configs](https://github.com/nf-core/configs).
 
---
 
## Citations
 
If you use WasteFlow 2.0 in your research, please cite:
 
> **WasteFlow 2.0**: [citation pending — DOI to be added on publication]
 
Please also cite the key underlying tools:
 
- **Nextflow**: Di Tommaso et al., *Nature Biotechnology* (2017). https://doi.org/10.1038/nbt.3820
- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **fastp**: Chen et al., *Bioinformatics* (2018). https://doi.org/10.1093/bioinformatics/bty560
- **Kraken2**: Wood et al., *Genome Biology* (2019). https://doi.org/10.1186/s13059-019-1891-0
- **Bracken**: Lu et al., *PeerJ Computer Science* (2017). https://doi.org/10.7717/peerj-cs.104
- **Krona**: Ondov et al., *BMC Bioinformatics* (2011). https://doi.org/10.1186/1471-2105-12-385
- **Minimap2**: Li, *Bioinformatics* (2018). https://doi.org/10.1093/bioinformatics/bty191
- **SAMtools**: Danecek et al., *GigaScience* (2021). https://doi.org/10.1093/gigascience/giab008
- **Mosdepth**: Pedersen & Quinlan, *Bioinformatics* (2018). https://doi.org/10.1093/bioinformatics/btx699
- **iVar**: Grubaugh et al., *Genome Biology* (2019). https://doi.org/10.1186/s13059-018-1618-7
- **snpEff**: Cingolani et al., *Fly* (2012). https://doi.org/10.4161/fly.19695
- **Freyja**: Karthikeyan et al., *Nature* (2022). https://doi.org/10.1038/s41586-022-05049-6
- **MultiQC**: Ewels et al., *Bioinformatics* (2016). https://doi.org/10.1093/bioinformatics/btw354
A full list of citations is available in [CITATIONS.md](CITATIONS.md).
 
---
 
## Contributing
 
Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to add new features, viral targets, or bug fixes.
 
For questions and support, please open a [GitHub Issue](https://github.com/your-org/wasteflow/issues).
 
---
 
## License
 
WasteFlow 2.0 is released under the [MIT License](LICENSE).
 