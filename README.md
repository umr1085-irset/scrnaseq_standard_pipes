# Preprocessing and processing scRNA-seq pipelines

This repository encapsulates two pipelines and their respective launch scripts. Please read the section on how to launch the scripts

## Preprocessing pipeline

* Input: Seurat object (.rds file), with sample information stored under `orig.ident` in the metadata
* Steps:
	* Load input Seurat file
	* Compute ribosomal and mitochondrial read percentages
	* Remove cells with less than 200 features
	* Remove genes with less than 10 expressing cells
	* Find cell outliers using `Scater's runColDataPCA()`, with ribosomal and mitochondrial read percentages, number of total reads and number of expressed genes as variables. Results stored under `scateroutlier`. PCA coordinates stored under `scateroutlierPC1` and `scateroutlierPC2`
	* Find doublet cells using `DoubletFinder`. This is performed on individual samples. Results stored under `DoubletFinder`
	* Apply `Seurat's SCTransform() and CellCycleScoring()` 
	* Compute dimensionality reductions (PCA, UMAP)
	* Find clusters, stored under `seurat_clusters`
	* Save Seurat object

## Processing pipeline
* Input: Preprocessed Seurat object (.rds file)
* Steps:
	* Load input preprocessed Seurat object
	* Remove outliers and doublets
	* Clean object (umap, pca, SCT assay, multiple metadata columns)
	* Remove genes with less than 10 expressing cells
	* Apply `Seurat's SCTransform() and CellCycleScoring()`
	* Compute dimensionality reductions (PCA, UMAP)
	* Find clusters
	* Save Seurat object

## Use launch scripts
The launch scripts are to be used with `Slurm's sbatch` command on a server. Each script has Slurm options that can modified:

```
#SBATCH --job-name="pre_pipe"
#SBATCH --output=pre_pipe.out
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --partition=sihp
#SBATCH --mail-user=<user-email>
#SBATCH --mail-type=ALL
#SBATCH --chdir=<output-dir>
```
* `job-name`: Name of Slurm job 
* `output`: Name of log file
* `mem`: Assigned memory (nust end with a `G` character)
* `cpus-per-task`: Number of assigned cps
* `partition`: cluster partition
* `mail-user`: user email for notification purposes
* `mail-type`: type of notifications (leave it set to `ALL`)
* `chdir`: Output directory where files will be saved

Users have to reference the pipeline input and output files directly in the launch scripts. The pipeline scripts should NOT be modified in that regard.

To launch a script:
```
sbatch pre_pipe_launcher.sh
```
