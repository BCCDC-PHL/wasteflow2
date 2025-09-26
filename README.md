# BCCDC-PHL/wasteflow2

## Introduction  
**BCCDC-PHL/wasteflow2** is a bioinformatics workflow designed for the processing and analysis of probe-enriched wastewater sequencing data for **viral surveillance**.  
The pipeline currently supports three viral targets of high public health relevance:  

- **SARS-CoV-2**  
- **Respiratory Syncytial Virus (RSV-A and RSV-B)**  
- **Influenza (A and B)**  

Built on [Nextflow DSL2](https://www.nextflow.io/) and the [nf-core](https://nf-co.re) framework, the workflow ensures scalability, reproducibility, and portability across different compute environments (local, HPC, cloud, Docker/Singularity).  

---

## Workflow Summary  

The pipeline consists of the following major steps:  

1. **Read Quality Control**  
   - Read trimming, adapter removal, and filtering using [`fastp`](https://github.com/OpenGene/fastp).  
   - Quality control reporting using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregated with [`MultiQC`](http://multiqc.info/).  

2. **Read Classification and Filtering**  
   - Reads are aligned against a **combined Viral + Human reference database**.  
   - Reads are classified and tagged with their respective **NCBI Taxonomy IDs**.  
   - **Human reads** are excluded from downstream analyses.  
   - Taxonomic profiles are generated using [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) and visualized interactively with [`Krona`](https://github.com/marbl/Krona/wiki).  
   - **Nextflow branching structure** is used to separate and filter reads per target taxon, enabling downstream parallel processing of SARS-CoV-2, RSV, and Influenza.  

3. **Per-Target Genome Processing**  

   ### SARS-CoV-2 and RSV (A & B)  
   - Reads mapped to their respective references using the `prepare_genome` subworkflow.  
   - RSV-A and RSV-B are first aligned jointly to an RSV reference database; BAM files are then **split into RSV-A and RSV-B** for independent processing.  
   - Each BAM undergoes:  
     - Indexing  
     - Coverage calculation and statistics  
   - Variant calling:  
     - Currently supported: [`iVar`](https://andersen-lab.github.io/ivar/html/)  
     - Planned: [`FreeBayes`](https://github.com/freebayes/freebayes)  
   - Variants from iVar are:  
     - Exported as `.tsv`, converted to VCF  
     - Annotated using [`SnpEff`](https://pcingola.github.io/SnpEff/) and [`SnpSift`](https://pcingola.github.io/SnpEff/se_snpSift/)  
     - Screened for **key mutations** based on a mutation watchlist  
   - VCFs are also processed with [`Freyja`](https://github.com/andersen-lab/Freyja) for lineage deconvolution.  

   ### Influenza  
   - (To be added — process differs slightly from the above workflow.)  

4. **Final Reporting**  
   - Results from all stages are collated and summarized into an integrated [`MultiQC`](http://multiqc.info/) report.  

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run BCCDC-PHL/wasteflow2 \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

BCCDC-PHL/wasteflow2 was originally written by Zohaib Anwar.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use BCCDC-PHL/wasteflow2 for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
