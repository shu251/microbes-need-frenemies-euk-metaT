# Microbes Need Frenemies
Descriptions and code for metatranscriptome analysis related to project "Microbes need frenemies"


## Eukrhythmic

[Instructions for running eukrhythmic](https://eukrhythmic.readthedocs.io/en/latest/index.html).

Running on HPC with slurm.

Followed all install instructions. See examples of set up documents, _config.yaml_ and _cluster.yaml_ for running snakemake workflow.

Use python submit scripts that come packaged with eukrhythmic. Using flag **--np** to perform dry run first.

```
conda activate eukrhythmic
python submit/eukrhythmic -np # to generate dry run, ensure everything runs ok.

# To run
python submit/eukrhythmic all
```


## Local & remote workflow for data analysis

### Make local changes

Add new changes to the git repo. Watch memory usage and keep large files separate. Add commit messages to annotate changes.
```
git add <newly-changed-files>
git commit -m "I changed my files!"
```

Return to remote:

Stay on branch on remote. But you can pull changes from main to your branch.
```
git pull origin main
```
Uploaded to a local staging area where you can render (For quarto github page) and keep modifying. When ready, push it to the main repository.
```
git push
```


### Make changes remotely

First make sure local repo and remote are up to date. 

```
git pull
```

You have the option to create a remote branch. Then merge the changes later.

```
git add <newly-changed-files>
git commit -m "I changed my files!"
```