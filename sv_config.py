#!/usr/bin/env python3
"""
For SV calling

"""
import os, sys
home = os.path.expanduser("~")
sys.path.append(f'{home}/jb/module')
sys.path.append(f'/pipeline')   # withiin singularity image
from sv_call import perform_sv_call

config = {
    'email':'huixiong9945@gmail.com',             # email address
    'pw':'/mnt/d/task/nsclc_wgs_batch1_4sample',                # working path
    'gm_build':'hg38',      # genome build, hg38 / hg19
    'sra_to_fastq' : 0,      # if start file is SRA file, set 1, if fastq file , set as 0
    'download_sra' : 0,      #if start point = SRR accession, but not yet downloaded, set as 1
    'paired' :1,            # read layout, if paired, set as 1

    # read group information
    # the file list could be fastq file, local SRA file or SRA accession ID
    # if SRA file or SRA accession, one file for each read_group :[ <file/accession>] pair
    'filelist' : {
'NCIH1435_LUNG_': ['/home/chenh19/tmp/fastq_wgs/batch1/SRR8652141/SRR8652141_1.fastq.gz',
                    '/home/chenh19/tmp/fastq_wgs/batch1/SRR8652141/SRR8652141_2.fastq.gz'],
 'NCIH1975_LUNG_': ['/home/chenh19/tmp/fastq_wgs/batch1/SRR8652050/SRR8652050_1.fastq.gz',
                    '/home/chenh19/tmp/fastq_wgs/batch1/SRR8652050/SRR8652050_2.fastq.gz'],
 'NCIH2030_LUNG_': ['/home/chenh19/tmp/fastq_wgs/batch1/SRR8652052/SRR8652052_1.fastq.gz',
                    '/home/chenh19/tmp/fastq_wgs/batch1/SRR8652052/SRR8652052_2.fastq.gz'],
 'NCIH522_LUNG_': ['/home/chenh19/tmp/fastq_wgs/batch1/SRR8670691/SRR8670691_1.fastq.gz',
                   '/home/chenh19/tmp/fastq_wgs/batch1/SRR8670691/SRR8670691_2.fastq.gz'],
        }
}

# the following is the default settings
config_default = {
    # 'dock_file': '/gpfs23/scratch/cqs/chenh19/dock/bioinfo.sif',    # if '', run the command from current env,
    'dock_file': '/home/chenh19/dock/bioinfo.sif',    # if '', run the command from current env,
    'read_group_from_fastq': 1,  # get the readgroup info from fastq, otherwise, would use the keys in "filelist" as the @RG
    'chunk': 100,  # split the fastq files into pieces, default = 50
    'no_cal': True,  # default = do recalibration after BWA( no_cal = False), if set to True, would skip the recal step
    'ignore_lb': ['call_sv', 'sra2fastq'], # skip summit script from certein steps, steps: ['sra2fastq ', 'split_fastq', 'bwa_align', 'bam_refine', 'bam_merge', 'call_sv']
    'region_bed_file': '', # extract the region from this bed files,  4th column would be the region name
    'tmp_dir': '/home/chenh19/tmp/fastq_wgs/batch1'   # tmp dir path for split fastq
}
perform_sv_call.run(config, config_default)

# copy this config to working path
this_script = os.path.abspath(__file__)
os.system(f'cp {this_script} {config["pw"]} 2>/dev/null')
