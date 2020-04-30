

## 0. intruduction

The scripts in this repository aim to automate the general SV call process.

The scripts are ready to be submitted to the ACCRE slurm system

There are 3 files.

### sv_config.py

The main content is 2 dict

config  -> the most frequently used parameters

config_default  -> the default settings, usually not need to be changed



here are the explanation for the keys

1. **email**

   email address, used for monitoring the status on ACCRE cluster

2. **pw**

   working path, all the results, scirpts and log would be stored under here, please make sure you have write permission here

3. **gm_build**

   genomic build,  currently only support hg19 / hg38

4. **sra_to_fastq**

   if initial file are SRA file, set 1 would run the fastq-dump process

   if initial  files are fastq format , set as 0

5. **download_sra**

   if start point = SRR accession(not downloaed to local disk), set as 1

6. **paired**

   if paired end, set as 1, otherwise 0

7. **filelist**

   key = label for the files. Also used as read group if read group is not available in fastq file

   value = fastq / SRA file / SRA accession, must be in a list even if there is only 1 file

8. **dock_file**

   Packing of essential softwares, if not available, just set as empty `''`

9. **read_group_from_fastq**

   if the read group is found in the fastq header, would extract this directly.

   If set as 0, would just use the key in filelist

10. **chunk**

    split the original fastq files into small chunks to accelerate the bwa mapping and reaglignment.

    default = 100
    
11. **no_cal**

    default = do recalibration after BWA( no_cal = False), 

    if set to True, would skip the recal step

12. **region_bed_file**

    extract the region from this bed files,  4th column would be the region name

    

### 

## 1. dependency

The softwares are packed inside of a singularity image file

```
/gpfs23/scratch/cqs/chenh19/dock/bioinfo.sif
```

You may skip installing the softwares below if you have access to the image above

### 1.1 softwares

1. python3.7  -> for run this pipeline

2. bwa  

3. samtools

4. gatk3

5. picard

6. sambamba  -> for sorting and merging

7. smoove  -> best practise of lumpy. 

   smoove requires python2.7

8. sra-toolkits

9. gnu parallel  

### 1.2 reference files

The reference files are available in ACCRE.

If you don't have access to the following files, you need to download these files first

and modify them in function `get_ref`  in `sv_call_modules.py`

```
# hg19
gm = /scratch/cqs_share/references/broad/hg19/v0/Homo_sapiens_assembly19.fasta
dbsnp = /scratch/cqs_share/references/broad/b37/dbsnp_138.b37.vcf.gz
mills_indel = /scratch/cqs_share/references/broad/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
excludebed=/gpfs23/scratch/h_vangard_1/chenh19/ref/svcall/ceph18.b37.lumpy.exclude.2014-01-15.bed


# hg38
gm = /scratch/cqs_share/references/broad/hg38/v0/Homo_sapiens_assembly38.fasta
dbsnp = /scratch/cqs_share/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
mills_indel = /scratch/cqs_share/references/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_indel = /scratch/cqs_share/references/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
excludebed = /gpfs23/scratch/h_vangard_1/chenh19/ref/svcall/exclude.cnvnator_100bp.GRCh38.20170403.bed
```



## 2.  example workflow

1. clone these repository

   ```
   cd /local/path
   git clone https://github.com/zeissmania/sv_tool.git
   ```

2. add the path above to PYTHONPATH 

   ```
   export PYTHONPATH=/local/path:$PYTHONPATH
   ```

3. edit the config file

4. run the config file, generate the pbs scripts

   ```
   python3 sv_call_config.py
   ```

5. open the working path set in config file, e.g. `/path/to/work`

   ```
   cd /path/to/work/final_submit
   
   # if submit to ACCRE
   bash final_submit.sh    
   
   
   # if run locally
   bash final_local_run.sh
   ```

   











