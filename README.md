# Linder et al 2022

This is a Nextflow pipeline for reproducing the analyses carried out in Linder et al, 2022. Code for the individual figures and tables is present in the `Figure_code` directory and are not yet a part of the nextflow pipeline. Code for recreating the data used in the figures and tables is present in the `bin` folder. The `helper_functions.R` file contains functions used in multiple steps of the pipeline and are sourced as needed. The minimal files needed to run this analysis are present in the `assets` folder, with the exception of the raw haplotype frequency file which can be found [here](https://wfitch.bio.uci.edu/~tdlong/sandvox/publications.html) in the `Linder et al 2022` link. This can be downloaded and copied into the `assets` folder to run the pipeline.  


[![homepage](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://nextflow.io/ "Redirect to nextflow homepage")

## Requirements

- Unix-like operating system (Linux, macOS, etc)
- Java 8 or later
- Docker 20.10.17 or later
- Nextflow 22.10.1 or later

## Quickstart

If you don't already have Nextflow, install by using the following command:

```
curl -s https://get.nextflow.io | bash
```

The Dockerfile and conda.yaml are included if updated or additional software tools are desired to be added to the Docker image used throughout this pipeline. To download the Docker image from dockerhub, use this command:

```
docker pull rlinder02/linder-et-al-2022:v1.1.0
```
An additional docker image from quay.io is also used for finding SNVs that potentially impact transcription factor binding sites that can be retrieved using this command:

```
docker pull quay.io/biocontainers/bioconductor-motifbreakr:2.12.0--r42hdfd78af_0
```

Clone the pipeline and then run on your local machine from GitHub:

```
git clone https://github.com/rlinder02/Linder-et-al-2022.git
```

```
nextflow run main.nf
```

To run the pipeline on AWS Batch, configure an AWS account to work with Nextflow as detailed [below](#pipeline-description) and then enter:

```
nextflow run main.nf -profile batch
```

The nextflow.config file ensures that the Docker image is used by default for containerised execution, unless overriden at the command line.


## Pipeline Description

This pipeline is used to analyze data presented in Linder et al, 2022. Currently, this pipeline is configured to be run either locally or on AWS. Input/output files can be stored either locally or in an S3 bucket. Instructions for configuring an AWS account to work with Nextflow can be found [here](https://staphb.org/resources/2020-04-29-nextflow_batch.html) and more generally [here](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-2-of-3/).

## Pipeline parameters

At the command line, you can specify several options (documented [here](https://www.nextflow.io/docs/latest/)). The most relevant options for this pipeline are listed as `params.` variables in the nextflow.config and main.nf file and can be specified on the command line line when launching Nextflow.

The number of CPUs and requested memory can be changed in the nextflow.config file and depends on the resources available. It is important to remember that, when running on AWS, the number of vCPUs provisioned for a job will not be the same as the number of CPUs reqested in the config file, as there are normally 2+ vCPUs per CPU (exact specifications for a variety of EC2 instance types are listed in table form [here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/cpu-options-supported-instances-values.html)).

## Pipeline results

The results are copied to a generic `results` folder, with subdirectories storing relevant files that can be used for further analyses and to examine the data in more detail.

## Credits

The general format of the pipeline was inspired by the work done [here](https://github.com/nextflow-io/rnaseq-nf) and [here](https://github.com/nf-core/rnaseq).



