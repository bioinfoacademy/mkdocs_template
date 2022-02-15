### A. Prepare samplesheet.csv file
1. Current [samplesheet.csv](../../Protocol/samplesheet.csv) file is a comma separated text file. This file could be edited in a plain text editor or a spreadsheet program. The file content is given below;
```
Lane,Sample,Index
*,WT_Gex,SI-GA-C1
*,Tg_Gex,SI-GA-D1
*,WT_Barcode,SI-GA-E1
*,Tg_Barcode,SI-GA-F1
```
In the Lane column, * specifies all lanes.
Sample column specifies sample names.
Index column specifies Index used for these samples (this info comes from [HS017_Run_info.pdf](../../Protocol/HS017_RUN_info.pdf))

### B. Generate fastq files from BCL files (AND Demultiplex at sample level)
1. Login to biowulf and request an interactive node using the below 'sinteractive' command, with 32 cores, 60gb RAM and 300gb of storage;
```
sinteractive --cpus-per-task=32 --mem=60g --gres=lscratch:300
```
2. Once the interacive node is assigned, run the below command to load the cellranger module;
```
# Change to the working directory
cd /data/Shih_Lab/scRNAseqNilisha
# Load the cellranger module
module load cellranger
```
3. After loading the cellranger module, run the below command to generate fastq files;
```
cellranger mkfastq --id=brainfastq --run=RawData/BCL --csv=Protocol/samplesheet.csv --localcores=$SLURM_CPUS_PER_TASK --localmem=34
```
The mkfastq is the subcommand to generate fastq files from BCL files.<br>
--id specifies the 'brainfastq' as a output folder (will be created new)<br>
--run specifies the location of the BCL folder that contains the raw BCL and associated files<br>
--csv specifies the samplesheet.csv file that we created in step A above<br>
--localcores specifies the number of cpus to use - same as that of the interactive node cpus<br>
--localmem specifies the 34 GB for memory

4. The mkfastq run takes about one hour to complete. Once the run is complete, the fastq files are available in the [outs/fastq_path/HCJY2DRXY/ folder within the output directory specified in the step #3](../../RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/).<br>

A preliminary Report of the flow cell summary and lane summary for each of the sample indexes is available in this [html](../../RawData/mkfastq-output/outs/fastq_path/Reports/html/index.html) file


### C. Cell count at individual mouse level
1. Prepare the multiplex_feature_reference.csv file
	* The multiplex_feature_reference.csv is a comman separated text file.
	* The multiplex_feature_reference.csv lists the hash tags (barcodes) used for each mouse.
	* Each mouse is labelled with a unique id (like WT_B1), a name (could be the same as id), read and pattern (which are standard for TotalSeqB protocol tags).
	* The multiplex_feature_rference.csv file for this run is located in the [Protocol](../../Protocol) folder. 
	* The contents of multiplex_feature_reference.csv for this run is given below.
```
id,name,read,pattern,sequence,feature_type
WT_B1,WT_B1,R2,5PNNNNNNNNNN(BC),ACCCACCAGTAAGAC,Multiplexing Capture
WT_B2,WT_B2,R2,5PNNNNNNNNNN(BC),GGTCGAGAGCATTCA,Multiplexing Capture
WT_B5,WT_B5,R2,5PNNNNNNNNNN(BC),CTTTGTCTTTGTGAG,Multiplexing Capture
WT_B6,WT_B6,R2,5PNNNNNNNNNN(BC),TATGCTGCCACGGTA,Multiplexing Capture
Tg_B3,Tg_B3,R2,5PNNNNNNNNNN(BC),CTTGCCGCATGTCAT,Multiplexing Capture
Tg_B4,Tg_B4,R2,5PNNNNNNNNNN(BC),AAAGCATTCTTCACG,Multiplexing Capture
Tg_B7,Tg_B7,R2,5PNNNNNNNNNN(BC),GAGTCTGCCAGTATC,Multiplexing Capture
Tg_B8,Tg_B8,R2,5PNNNNNNNNNN(BC),TATAGAACGCCAGGC,Multiplexing Capture
```
2. Prepare the multiplex_config_new.csv file
	* The multiplex_config.csv is text file, with three sections (gene-expression, libraries, samples).
	* The gene-expression section has information about the reference genome location, expected number of cells and the hash-tag (cmo-set) reference file (as described in #1 above) location
	* The libraries section describes the fastq files folders for each of the sample
	* The sample section describes the sample ids for each of the mouse (can use the same ids for sample_ids, cmo_ids and description
	* The file used for this analysis is available in the [Protocol](../../Protocol) folder
	* The contents of the multiplex_config.csv file used for this analysis is given below
```
[gene-expression]
reference,/fdb/cellranger/refdata-gex-mm10-2020-A
expect-cells,6000
cmo-set,/data/Shih_Lab/scRNAseqNilisha/Protocol/multiplex_feature_reference.csv

[libraries]
fastq_id,fastqs,feature_types
WT_Gex,/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/WT_Gex,Gene Expression
Tg_Gex,/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/Tg_Gex,Gene Expression
WT_Barcode,/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/WT_Barcode,Multiplexing Capture
Tg_Barcode,/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/Tg_Barcode,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
WT_B1,WT_B1,WT_B1
WT_B2,WT_B2,WT_B2
WT_B5,WT_B5,WT_B5
WT_B6,WT_B6,WT_B6
Tg_B3,Tg_B3,Tg_B3
Tg_B4,Tg_B4,Tg_B4
Tg_B7,Tg_B7,Tg_B7
Tg_B8,Tg_B8,Tg_B8
```
3. Login to biowulf and request an interactive node using the below 'sinteractive' command, with 32 cores, 64gb RAM and 300gb of storage;
```
sinteractive --cpus-per-task=32 --mem=64g --gres=lscratch:300
```
4. Once the interacive node is assigned, run the below command to load the cellranger module;
```
# Change to the working directory
cd /data/Shih_Lab/scRNAseqNilisha
# Load the cellranger module
module load cellranger
```
5. After loading the cellranger module, run the below command to generate cell counts;
```
cellranger multi --id=wttgbrain_hto_demultiplexed --csv=Protocol/multiplex_config.csv --jobmode=slurm --maxjobs=25
```
The multi is the subcommand to demultiplex based on HTO samples.<br>
--id specifies the 'wttgbrain_hto_demultiplexed' as a output folder (will be created new)<br>
--csv specifies the location of the multiplex_config.csv file (Protocol/multiplex_config.csv)<br>
--jobmode=slurm specifies that the job needs to be run in batch mode<br>
--maxjobs specifies the number of jobs allowed for this batch run (25 in our run)<br>

6. The cellranger multi runs for a few hours to complete. Once the run is complete, the output files are available in the [outs/per_sample_outs folder within the 'wttgbrain_hto_demultiplexed' directory specified in the step #5](../../Results/wttgbrain_hto_demultiplexed/outs/per_sample_outs).<br>

A preliminary Report with the count summary, alignment quality, some preliminary analysis are available in the [per_sample_outs](../../Results/wttgbrain_hto_demultiplexed/outs/per_sample_outs) folder, within each of the sample folders. For example, this [web_summary.html](../../Results/wttgbrain_hto_demultiplexed/outs/per_sample_outs/Tg_B3/web_summary.html) file is available in the [per_sample_outs/Tg_B3](../../Results/wttgbrain_hto_demultiplexed/outs/per_sample_outs/Tg_B3) folder<br>
The matrix files that we will use in the downstream analysis in seurat is also available in the per samples folder.

### D. Cell count at sample level
1. Prepare libraryTG.csv and libraryWT.csv files (these are comma separated text files)
	* The libraryTG.csv describes the location of the fastq files folder for the TG sample, the corresponding sample types (Tg) and the library preparation type (gene expression or antibody capture)
	* The libraryTG.csv for this run is located [here](../../Protocol/libraryTG.csv) in the Protocol [folder](../../Protocol)
	* Similarly, prepare the libraryWT.csv file for the WT sample.
	* The contents of the libraryTG.csv file and libraryWT.csv are given below. 

!!! note
	Since we are going to the data only at the samples level (wt vs tg), we do not have to include the antibody capture information in this run, but we are including it just in case.
```
fastqs,sample,library_type
/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/WT_Barcode/,WT_Barcode,Antibody Capture
/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/WT_Gex/,WT_Gex,Gene Expression
```

```
fastqs,sample,library_type
/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/Tg_Barcode/,Tg_Barcode,Antibody Capture
/data/Shih_Lab/scRNAseqNilisha/RawData/mkfastq-output/outs/fastq_path/HCJY2DRXY/Tg_Gex/,Tg_Gex,Gene Expression
```
2. Prepare a common feature_reference.csv file(this is a comma separated text file)<br>

* The feature_reference.csv file describes the antibodies (hash tag) used to hash tag the individual mouse (like B1, B2 etc.), a name for the individual mouse (like WT_B1), and then the pattern R2,5PNNNNNNNNNN(BC) (this is specific to [TotalSeq B protocol](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis), the actual Tag barcode sequence and the feature type (antibody capture)<br>
* The feature_reference.csv for this run is located [here](../../Protocol) in the Protocol [folder](../../Protocol)<br>
* The contents of the feature_reference.csv file is given below.<br>

```
id,name,read,pattern,sequence,feature_type
WT_B1,WT_B1,R2,5PNNNNNNNNNN(BC),ACCCACCAGTAAGAC,Antibody Capture
WT_B2,WT_B2,R2,5PNNNNNNNNNN(BC),GGTCGAGAGCATTCA,Antibody Capture
WT_B5,WT_B5,R2,5PNNNNNNNNNN(BC),CTTTGTCTTTGTGAG,Antibody Capture
WT_B6,WT_B6,R2,5PNNNNNNNNNN(BC),TATGCTGCCACGGTA,Antibody Capture
Tg_B3,Tg_B3,R2,5PNNNNNNNNNN(BC),CTTGCCGCATGTCAT,Antibody Capture
Tg_B4,Tg_B4,R2,5PNNNNNNNNNN(BC),AAAGCATTCTTCACG,Antibody Capture
Tg_B7,Tg_B7,R2,5PNNNNNNNNNN(BC),GAGTCTGCCAGTATC,Antibody Capture
Tg_B8,Tg_B8,R2,5PNNNNNNNNNN(BC),TATAGAACGCCAGGC,Antibody Capture
```
3. Prepare the swarm file containing the commands for both the samples (WT and TG)<br>

* This swarm file (saved as [count.swarm](../../Tools)) is a text file with the below commands.
* The swarm file specifies the 'count' subcommand for cellranger, the output folder name (--id), path to the library file (--libraries), path to the reference transcriptome (--transcriptome), path to the feature reference file (--feature-ref), the expected number of cells (--expect-cells), the number of cores to use and the amount of memory to use.

```
cellranger count --id=WTbrain --libraries=/data/Shih_Lab/scRNAseqNilisha/Protocol/libraryWT.csv --transcriptome=/fdb/cellranger/refdata-gex-mm10-2020-A --feature-ref=/data/Shih_Lab/scRNAseqNilisha/Protocol/feature_reference_WT.csv --expect-cells=6000 --localcores=28 --localmem=64

cellranger count --id=TGbrain --libraries=/data/Shih_Lab/scRNAseqNilisha/Protocol/libraryTG.csv --transcriptome=/fdb/cellranger/refdata-gex-mm10-2020-A --feature-ref=/data/Shih_Lab/scRNAseqNilisha/Protocol/feature_reference_TG.csv --expect-cells=6000 --localcores=28 --localmem=64
```
4. Change to the Tools directory using the below command<br>

```
cd /data/Shih_Lab/scRNAseqNilisha/Tools/
```
5. Run the swarm job using the below command, specifying the swarm file name, memory (-g), threads (-t), scratch space (--gres) and the module cellranger;

```
swarm -f count.swarm -g 64 -t 28 --gres=lscratch:200 --module cellranger
```

!!! note
    Note down your job id.

Check your job status using the "jobload" command or the squeue command;

!!! note
    Replace the job id in the below command with your job id

```
squeue --job 14354599
```
6. Extend the swarm job running time using the below command;

!!! note
    Replace the job id in the below command with your job id

```
newwall --jobid 13558558 --time 36:00:00
```
7. Output<br>

* The cellranger count runs for a few hours to complete. Once the run is complete, the output files are available in the outs folder within the 'WTbrain' and 'TGbrain' directories specified in the swarm file. In this run we moved the WTbrain and TGbrain results folders to the [Results folder](../../Results).<br>
* A preliminary Report with the count summary, alignment quality, some preliminary analysis are available in the respective outs folder in the [web_summary.html file](../../Results/WTbrain/outs/web_summary.html)<br>
* The respective clustering file [(cloupe.cloupe)](../../Results/WTbrain/outs/) based on the preliminary analysis is also available in the 'outs' folder, to visualize and explore the t-SNE and UMAP data.<br>
* The matrix files that we will use in the downstream analysis in seurat is also available in the 'outs' folder.

### E. Combine Cell counts (using sample level counts), to get a combined loupe file
1. Prepare the aggr.csv file

	* The aggr.csv is a comma separated text file.
	* The aggr.csv lists the sample Id, path to the molecule_h5 file for each of the samples obtained from the Step D above (WT and TG) and the sample group (WT or TG)
	* The aggr.csv file for this run is located in the [Tools](../../Tools) folder. 
	* The contents of aggr.csv for this run is given below.
```
sample_id,molecule_h5,SampleGroup
WT,/data/Shih_Lab/scRNAseqNilisha/Results/WTbrain/outs/molecule_info.h5,WT
TG,/data/Shih_Lab/scRNAseqNilisha/Results/TGbrain/outs/molecule_info.h5,TG
```

2. Login to biowulf and request an interactive node using the below 'sinteractive' command, with 32 cores, 64gb RAM and 300gb of storage;
```
sinteractive --cpus-per-task=32 --mem=64g --gres=lscratch:300
```
3. Once the interacive node is assigned, run the below command to load the cellranger module;
```
# Change to the working directory
cd /data/Shih_Lab/scRNAseqNilisha
# Load the cellranger module
module load cellranger
```
4. After loading the cellranger module, run the below command to combine the cell counts from WT and TG and generate the combined loupe file;
```
cellranger aggr --id=WT_TG_combined --csv=Tools/aggr.csv
```
The aggr is the subcommand to combine the sample counts from WT and TG results folders.<br>
--id specifies the 'WT_TG' as a output folder (will be created new)<br>
--csv specifies the location of the aggr.csv file (Tools/aggr.csv)<br>

5. The cellranger aggr runs for a few hours to complete. Once the run is complete, the output files are available in the [outs folder within the 'WT_TG_combined' directory specified in the step #4 above](../../Results/WT_TG_combined/outs/).<br>

A preliminary Report with the count summary, alignment quality, some preliminary analysis are available in the [web_summary.html](../../Results/WT_TG_combined/outs/web_summary.html) file is available in the [outs](../../Results/WT_TG_combined/outs) folder<br>
The combined loupe file is available in the [count](../../Results/WT_TG_combined/outs/count) folder

### Software
1. cellranger  6.0.0
2. R 4.0.5
3. mkdocs 1.1.2

### References

1. [fastq generation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)
2. [Demultiplex and count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi)
2. [TotalSeq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct)
3. [CellRanger Biowulf](https://hpc.nih.gov/apps/cellranger.html)
4. [Swarm Biowulf](https://hpc.nih.gov/apps/swarm.html)




