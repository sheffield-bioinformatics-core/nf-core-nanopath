# ![nf-core/nanopath](docs/images/nf-core-nanopath_logo_light.png#gh-light-mode-only) ![nf-core/nanopath](docs/images/nf-core-nanopath_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanopath/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/nanopath)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nanopath-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nanopath)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/nanopath** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Initialize the data:
      If a fastq directory is provided:
         Concatenate fastq files using CAT_FASTQS.

2. Validate input:
      Use the INPUT_CHECK subworkflow to read samplesheet, validate, and stage input files.
      Branch reads based on their status (discontinued or samples).

3. Perform Quality Control:
      Run ([`FASTP`](https://github.com/OpenGene/fastp)) for quality control, filtering, and preprocessing.
      Filter out samples with no reads left after FASTP.
      Run ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) on the processed reads.

4. Classfy and Cluster:
      If specified, remove unclassified reads using ([`KRAKEN2`](https://github.com/DerrickWood/kraken2)).
      Subset reads based on specified parameters (default 100k reads to keep memory requirements reasonable).
      Perform k-mer frequency analysis with KMER_FREQS.
      Perform read clustering with READ_CLUSTERING using ([`HDBSCAN`](https://github.com/scikit-learn-contrib/hdbscan)) and ([`UMAP`](https://umap-learn.readthedocs.io/en/latest/)).

5. Split Clusters and Correct Errors:
      Split clusters.
      Perform error correction using ([`CANU`](https://github.com/marbl/canu)).

6. Select and Polish Draft:
      Select draft reads using ([`FASTANI`](https://github.com/ParBLiSS/FastANI)).
      Polish drafts using ([`RACON`](https://github.com/isovic/racon)).
      Generate final consensus using ([`MEDAKA`](https://github.com/nanoporetech/medaka)).

7. Classify Taxonomically:
      Based on chosen tool, classify consensus sequences with ([`BLAST`](https://www.ncbi.nlm.nih.gov/books/NBK279690/)), ([`SEQMATCH`](https://github.com/rdpstaff/SequenceMatch)), ([`KRAKEN`](https://github.com/DerrickWood/kraken2)) or all of them. 
      Join classification results using JOIN_RESULTS.

8. Estimate Abundace:
      Estimate abundance per sample per detected species. 

9. Generate Reports:
      If report generation is chosen:
         Generate HTML reports.
## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

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
nextflow run nf-core/nanopath \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://nf-co.re/nanopath/usage) and the [parameter documentation](https://nf-co.re/nanopath/parameters).

## Pipeline output

To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/nanopath/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/nanopath/output).

## Credits

nf-core/nanopath was originally written by Magdalena Dabrowska.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nanopath` channel](https://nfcore.slack.com/channels/nanopath) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/nanopath for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
