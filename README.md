# Synopsis

This repository contains the source code and input data used for the data analysis and modeling for XXX <!-- TODO: reference paper -->.

Input data download and subsequent analyses are automated using Nextflow and Singularity/Apptainer.

## Docker images availability

Docker images are hosted on GHCR and can be found in the [Packages](https://github.com/orgs/paulsengroup/packages?repo_name=2022-mcf10a-cancer-progression) page of this repository.

Images were generated using the `build*dockerfile.yml` GHA workflows using the Dockerfiles from the `containers` folder.

## Nextflow workflows

Nextflow workflows under `workflows` were developed and tested using Nextflow v22.10.7, and should in principle work with any version supporting Nextflow DSL2.

Each workflow is paired with a config file (see `configs` folder). As an example, `workflows/fetch_data.nf` is paired with config `configs/fetch_data.config`.

## Requirements

- Access to an internet connection (required to download input files and Docker images)
- Nextflow v20.07.1 or newer
- Apptainer/Singularity (tested with Apptainer v1.1.6)

## Running workflows

Please make sure Nextflow is properly installed and configured before running any of the workflows.

The following workflows should be executed first, as they download and prepare files required by other workflows.
1. `fetch_data.nf`
2. `preprocess_data.nf`

The `fetch_data.nf` workflow requires internet access and can fail for various reason (e.g. connection reset by peer, service temporarily unavailable etc.). In case the workflow fails, wait few minutes, then relaunch the workflow.

The execution order of the rest of the worklows varies depending on which parts of the data analysis you are interested in re-running.
The following order assumes you want to re-run the entire analysis. If you only want to re-run some steps, feel free to get in touch with us to know which steps you have to run.

1. `run_nfcore_hic.sh`
2. `run_postprocess_nfcore_hic.sh`
3. `run_tad_analysis.sh`
4. `run_compartment_analysis.sh`
5. `run_call_tad_cliques.sh`
6. `run_chrom3d.sh`
7. `run_nfcore_rnaseq.sh`
8. `run_diff_expression_analysis.sh`
8. `run_nfcore_chipseq.sh`
9. `run_nfcore_nascent.sh`
10. `run_comparative_analysis_hic.sh`
11. `run_comparative_analysis.sh`

Inside the `config` folder there are the following base configs:
- `base_hovig.config`
- `base_linux.config`
- `base_saga.config`
- `base_macos.config`

These configs are passed to all workflows, and define available computation resources.
You will most likely have to update one of the configs with resources available on our machine/cluster.
