# Microbes Need Frenemies
Descriptions and code for metatranscriptome analysis related to project "Microbes need frenemies"


## Eukrhythmic

[Instructions for running eukrhythmic](https://eukrhythmic.readthedocs.io/en/latest/index.html).

Running on HPC with slurm.

Followed all install instructions. See examples of set up documents, _config.yaml_ and _cluster.yaml_ for running snakemake workflow.

Use python submit scripts that come packaged with eukrhythmic. Using flag **--np** to perform dry run first.

```
python submit/eukrhythmic -np
```