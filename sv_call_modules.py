"""
the modules for sv_call pipeline
"""
import os
import re
import time

import urllib.request
import xml.etree.ElementTree as Et
from time import time as t
import logging

logger = logging.getLogger('main')

def validate_fastq(data):
    """
    input = 4 line of a candidate fastq files, type=str
    """
    data = newlist(data, '\n')
    if len(data) != 4:
        return 0
    try:
        line1 = 1 if data[0][0] == '@' else 0
        line2 = 1 if re.match('^[ATCGN-]+$', data[1]) else 0
        line3 = 1 if data[2][0] == '+' else 0
    except:
        return 0
    if all([line1, line2, line3]):
        return 1
    return 0


def validate_sv_config(pw, gm_build, sra_convert, filelist, download_sra, dock_file, ref_info, config_default):
    """
    validate the parameters

    """
    logger.info('validate sv config')
    err = []
    if gm_build not in ['hg19', 'hg38', 'grch37', 'grch38']:
        msg = f'wrong genome build, valid should be hg19 / hg38'
        logger.fatal(msg)
        err.append(msg)

    mds(f'{pw}/validate')

    # check the dock_file
    if dock_file:
        if not os.path.exists(dock_file):
            msg = f'dock_file: < {dock_file}  > not exist'
            logger.fatal(msg)
            err.append(msg)
    try:
        gm = ref_info['gm']
        mills = ref_info['mills_indel']
        dbsnp = ref_info['dbsnp']
        known_indel = ref_info['known_indel']
        if not os.path.exists(gm):
            msg = f'genome fastq file not found!\nthe expected path={gm}'
            logger.warning(msg)
            print('*'*20)
            gm_new = input(f'specify the genome fasta file for {gm_build}\n')
            if os.path.exists(gm_new):
                isize = os.path.getsize(gm_new)
                isize = isize / 1024 / 1024 / 1024
                logger.debug(f'new input gm fa({gm_new}): size = {isize}')
                if isize < 2:
                    logger.fatal(f'new input fasta file found, but filesize too small: {isize}G')
                    error.append('gm fasta file not found')
                else:
                    ref_info['gm'] = gm_new
                    logger.info(f'your input is {gm_new}, seems ok')
            else:
                logger.fatal(f'your new input also does not exist: {gm_new}')
                err.append(msg)

        for i in [mills, dbsnp, known_indel]:
            i = i.strip()
            if not os.path.exists(i) and not config_default['no_cal']:
                logger.warning('file not exist: {i}')
                err.append(f'file not exist: {i}')
    except:
        msg = f'key in the ref dict is wrong {ref_info}'
        logger.fatal(msg)
        err.append(msg)


    for lb, fls in filelist.items():
        for ifl in fls:
            if not os.path.isfile(ifl) and not download_sra:
                err.append(f'file not exist- {lb} : {ifl}')
            elif os.path.isfile(ifl):
                tmp = getlb(ifl)
                if sra_convert and not os.path.isfile(f'{pw}/validate/validated_{tmp}'):
                    fd_res = os.popen(f'fastq-dump -Z -X 1 {ifl} 2>/dev/null|wc -l').read()
                    if fd_res == '0':
                        err.append(f'wrong SRA format - {lb} : {ifl} ')
                    elif fd_res == '4':
                        os.system(f'touch {pw}/validate/validated_{tmp}')
                elif not sra_convert and not os.path.isfile(f'{pw}/validate/validated_{tmp}'):
                    ext = ifl.rsplit('.', 1)[-1]
                    if ext == 'gz':
                        data = os.popen(f'gunzip -c {ifl} 2>/dev/null|head -4').read().strip()
                        fastq_res = validate_fastq(data)
                    elif ext in ['fastq', 'fq']:
                        data = os.popen(f'head -4 {ifl} 2>/dev/null').read().strip()
                        fastq_res = validate_fastq(data)
                    else:
                        fastq_res = 0
                        err.append(f'unknown fastq extension, fastq/fq/gz - {lb} : {ifl}')
                    if fastq_res:
                        os.system(f'touch {pw}/validate/validated_{tmp}')
                    else:
                        err.append(f'wrong fastq format - {lb} : {ifl}')

    return err if err else 0


def sra_ebi_query(sra_acc, output):
    part1 = sra_acc[:6]
    part2 = f'00{sra_acc[-1]}'
    ftp_folder = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{part1}/{part2}/{sra_acc}'
    os.system(
        f'/home/chenh19/miniconda/bin/lftp ftp://ftp.sra.ebi.ac.uk <<< "cd {ftp_folder}/;nlist" >{output} 2>&1')
    flist = open(output).read().split('\n')
    return flist


def sra2fastq(pw, lb, sra_file, paired, email, dock_file):
    """
    pw = str, output fastq path
    sra_file = str, sra_file or just the SRA-accession
    paired = bool, sequence layout
    """

    s = t()
    time_count = {f'pos{_}': 0 for _ in range(1, 5)}
    time_count['ftp'] = 0

    step = 'sra2fastq'
    mds(f'{pw}/{step}', 'pbs,log,result,check')

    # fastq-dump commands
    split = '--split-e' if paired else ''

    m = re.match(r'.*([SED]RR\d+)', sra_file)
    assert m, print('wrong sra file name (SRR, ERR, DRR)', sra_file)

    time_count['pos1'] += t() - s
    s = t()

    sra_acc = m.group(1)

    # first, we try to find if the fastq file already exist in ebi server
    do_fastq_dump = 0

    ftp_query_result = f'{pw}/{step}/ftp_query_{lb}'

    part1 = sra_acc[:6]
    part2 = f'00{sra_acc[-1]}'
    ftp_folder = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{part1}/{part2}/{sra_acc}'

    try:
        if os.path.getsize(ftp_query_result) > 5:  # this file is solid
            flist = open(ftp_query_result).read().split('\n')
        else:
            flist = sra_ebi_query(sra_acc, ftp_query_result)
    except:
        flist = sra_ebi_query(sra_acc, ftp_query_result)
    flist = [_.strip() for _ in flist if _.strip()]

    time_count['ftp'] += t() - s
    s = t()

    if flist[0].find('Access failed') > -1:
        print('no fastq found on EBI server: would skip this sample', sra_file, lb)
        return '1'
        # do_fastq_dump = 1
    elif len(flist) not in [1, 2]:
        print(f'wrong fastq count on EBI n = {len(flist)}, would skip this sample', ' '.join(flist))
        return '2'
        # do_fastq_dump = 1
    else:
        outpw = f'{pw}/{step}/result/{sra_acc}'
        mds(outpw)
        url = ' '.join([f'{ftp_folder}/{ifl}' for ifl in flist])
        cmd = f"err=0;cd {outpw}\nparallel 'wget -c {{}}' ::: {url};"
        for ifl in flist:
            cmd += f"""

# check the fastq integrity
size_local=$(ls {outpw}/{ifl} -l|awk '{{print $5}}')
size_ebi=$(/home/chenh19/miniconda/bin/lftp ftp://ftp.sra.ebi.ac.uk <<< "cd {ftp_folder}/;ls {ifl}"|awk '{{print $5}}')
if [ $size_local -ne $size_ebi ];then
    echo incomplete fastq from EBI, {ifl}, size_local=$size_local,size_ebi=$size_ebi >> {pw}/{step}/check/{lb}.{step}.err
    err=1
fi
        """
        result = [f'{outpw}/{ifl}' for ifl in flist]
        # print('will download sra fastq files from EBI, expected file = ', result)

    if do_fastq_dump:
        # result list
        result = [f'{pw}/{step}/result/{sra_acc}_1.fastq.gz',
                  f'{pw}/{step}/result/{sra_acc}_2.fastq.gz'] if paired else [f'{pw}/{step}/result/{sra_acc}.fastq.gz']
        cmd = f'fastq-dump --gzip --helicos {split} --origfmt -O {pw}/{step}/result {sra_file}'

    time_count['pos2'] += t() - s
    s = t()

    # build pbs.sh
    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash
{cmd}

date

if [ -s {result[0]} ] && [ $err -eq 0 ];then
    echo 'completed@@@@@@@'
    touch {pw}/{step}/check/{lb}.{step}.done
fi
        """, file=out)
    # build pbs

    fl_pbs = build_pbs(pw, step, result[0], lb, email, dock_file, mem='40G')

    time_count['pos3'] += t() - s
    s = t()

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done', time_count]


def get_total_line(fastq_file_name, lb):
    m = re.match(r'.+_?(SRR\d+)_', fastq_file_name)
    if m:
        srr = m.group(1)
        try:
            reads = int(get_reads_from_ncbi(srr))   # if can't get the reads, return err, this section would fail
            totalline = reads * 4
            # print(f'get reads count from NCBI, label={lb}, SRR={srr}, reads={reads}')
        except:
            totalline = f'$(gunzip -c {fastq_file_name}|wc -l)'
    else:
        totalline = f'$(gunzip -c {fastq_file_name}|wc -l)'
    return totalline


def split_fastq(pw, lb, fastq_files, email, dock_file, chunk=50, tmp_root='', compress=False):
    """
    prepare for the bwa step
    the bwa mapping is very time consuming, spliting the fastq would accelerate the process
    these files would be written to the /tmp
    chunk = pieces number for analysis
    lb = label for the output

    """
     # extra split parameter if need to compress the file here, default is compress at the bwa step
    gzip_split_file = r""" --filter "gzip >\$FILE.done.gz " """ if compress else ''
    gz_flag = '.gz' if compress else ''
    tmp_root = tmp_root if tmp_root else f'{pw}/tmp'
    step = 'split_fastq'
    pwtmp = f'{tmp_root}/{step}'
    mds(pwtmp)
    mds(f'{pw}/{step}', 'pbs,log,check')

    # the estimate reads count step is extreamly slow if the fastq is from WGS
    # so we just record the reads count when the count is already done
    # one way is just wc -l the fastq file, but this is extremely slow
    # we can query the NCBI database, and get the reads count for this SRA,

    line_count_record = f'{pw}/{step}/line_count_{lb}'
    try:
        totalline = int(open(line_count_record).read().split('\n')[0].strip())
        if totalline < 100000:
            totalline = get_total_line(fastq_files[0], lb)
    except:
        totalline = get_total_line(fastq_files[0], lb)

    assert totalline > 100000, 'FastQtotal line less than 100k, please check'
    with open(line_count_record, 'w') as out:
        print(totalline, file=out)

    # result list
    result = []
    fq1 = fastq_files[0]
    fn_short = getlb(fq1)
    result1 = [f'{pwtmp}/{fn_short}_chunk/x{str(i).rjust(3, "0")}.done{gz_flag}' for i in range(chunk)]

    if len(fastq_files) == 2:
        fq2 = fastq_files[1]
        fn_short = getlb(fq2)
        result2 = [f'{pwtmp}/{fn_short}_chunk/x{str(i).rjust(3, "0")}.done{gz_flag}' for i in range(chunk)]
        result = list(zip(result1, result2))
    else:
        result = [[_] for _ in result1]

    # build pbs.sh
    # print(result[0][-1])
    # determine read2
    if len(fastq_files) == 2:
        read2 = fastq_files[1]
    elif len(fastq_files) == 1:
        read2 = ''

    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash

# get the total line of this fastq file
cd {pwtmp}
totalline={totalline}
echo $totalline >{line_count_record}
# reads_per_file=$(echo "$totalline/{chunk}/4 +1"|bc)
reads_per_file=$(perl -e "print int($totalline/{chunk}/4) +1")
# lines_per_file=$(echo "$reads_per_file*4"|bc)
lines_per_file=$(perl -e "print $reads_per_file*4")
export lines_per_file
read1={fastq_files[0]}
read2={read2}

""", file=out)
        donefile = f"{pw}/{step}/check/{lb}.{step}.done"
        print(
            """

# if the final filechunk is not %s, just abort, would not the gzip and rename step
# here we use parallel to run splitting of fastq files from 2 reads simultaneously

# do gzip in following step
parallel --env lines_per_file 'x={/};x=${x%%.gz};\\
    x=${x%%.fastq};x=${x%%.fq};\\
    if [ -d ${x}_chunk ];then rm -rf ${x}_chunk;fi;\\
    mkdir ${x}_chunk;\\
    cd ${x}_chunk;\\
    split -d -a 3 -l $lines_per_file  %s <(gunzip -c {}) && echo "split is done";\\
    if [ `ls x*|wc -l` -ne %s ];then \\
        echo "wrong split chunk number:{}"; \\
        rm x*;exit 1;\\
    else \\
        ls x*|egrep -v "*.(gz|done)$"|xargs -I ,,, mv ,,, ,,,.done;\\
    fi' \\
::: $read1 $read2

date

if [ -s %s ];then
    echo 'completed@@@@@@@'
    touch %s
fi

""" % (chunk, gzip_split_file, chunk, result[0][-1], donefile), file=out)

    # build pbs
    fl_pbs = build_pbs(pw, step, result1[-1], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done']


def bwa_align(pw, split_raw, genome, email, read_group, dock_file, tmp_root=''):
    """
    if run the bwa from dock, just add the dock prefix, such as singularity exec <doc_img>
    """

    step = "bwa_align"
    tmp_root = tmp_root if tmp_root else f'{pw}/tmp'
    pwtmp = f'{tmp_root}/sambamba_sort'
    mds(pwtmp)
    mds(f'{pw}/{step}', 'pbs,log,result,check')

    fastq_files = []
    for i in split_raw:
        a, fn_short, split_i = i.rsplit('/', 2)
        split_i = split_i.split('.')[0]
        fastq_files.append('/'.join([a, fn_short, fn_short + split_i]) + '.gz')

    lb = getlb(fastq_files[0], rmsn=1)

    # get result
    result = [f'{pw}/{step}/result/{lb}.sorted.bam']

    # build pbs.sh
    mv2 = f'ln -s {split_raw[1]}.gz {fastq_files[1]} 2>/dev/null' if len(split_raw) == 2 else ''
    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash

if [ ! -s "{split_raw[0]}.gz" ];then
    gzip {" ".join(split_raw)} 2>/dev/null
fi
ln -s {split_raw[0]}.gz {fastq_files[0]} 2>/dev/null
{mv2}
echo bwa mem start
date
bwa mem -M -R "{read_group}" \\
    {genome} \\
    {" ".join(fastq_files)}|\\
    samtools view -bS |\\
    sambamba sort -m 20G \\
        --tmpdir {pwtmp} \\
        -o {pw}/{step}/result/{lb}.sorted.bam \\
        /dev/stdin
date

if [ -s {result[0]} ];then
    echo 'completed@@@@@@@'
    touch {pw}/{step}/check/{lb}.{step}.done
fi
        """, file=out)

    # build pbs
    fl_pbs = build_pbs(pw, step, result[0], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done']


def refine_bam(pw, bam_in, ref_info, email, dock_file, tmpforrg, tmp_root=''):
    """
    1. use gatk IndelRealigner to reduce the miscall of indels
    2. use BQSR to  reduce the effects of analysis artefacts produced by your sequencing machines
    3. mark dup
    """
    step = 'bam_refine'
    tmp_root = tmp_root if tmp_root else f'{pw}/tmp'
    pwtmp = f'{tmp_root}/gatk'
    mds(pwtmp)
    mds(f'{pw}/{step}', 'pbs,log,result,check')
    lb = getlb(bam_in, 1)

    genome = ref_info['gm']  # gm, dbsnp, mills_indel, known_indel
    mills = ref_info['mills_indel']
    dbsnp = ref_info['dbsnp']
    known_indel = ref_info['known_indel']
    known_sites_indel = f'-knownSites {known_indel}' if known_indel else ''
    known_indel = f'-known {known_indel}' if known_indel else ''

    # filename
    f_recal_indel_interval = f'{pw}/{step}/result/{lb}.indel.intervals'
    f_recal_snp_table = f'{pw}/{step}/result/{lb}.BQSR.table'
    bam_recal_indel = f'{pwtmp}/{lb}.recal.indel.bam'
    bam_add_rg = f'{pwtmp}/{lb}.addrg.bam'
    bam_recal_snp = f'{pwtmp}/{lb}.recal.dbsnp.bam'
    bam_markdup = f'{pw}/{step}/result/{lb}.recal.markdup.bam'
    filedone = f'{pw}/{step}/check/{lb}.{step}.done'

    # get result
    result = [bam_markdup]

    # build pbs.sh
    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash
# realign indel
echo "gatk RealignerTargetCreator"
date
if [ ! -s {f_recal_indel_interval} ];then
    gatk3 -Xmx2g -T RealignerTargetCreator \\
        -R {genome} \\
        -I {bam_in} \\
        -known {mills} {known_indel} \\
        -o {f_recal_indel_interval}
else
    echo {lb}.indel.intervals already exist {f_recal_indel_interval}
fi

echo -e "\\n\\ngatk IndelRealigner"
date
if [ ! -s {bam_recal_indel} ];then
    gatk3 -Xmx4g -T IndelRealigner \\
        -R {genome} \\
        -I {bam_in} \\
        -targetIntervals {f_recal_indel_interval} \\
        -known {mills} {known_indel} \\
        -o {bam_recal_indel}
else
    echo file already exist {bam_recal_indel}
fi

# temp, add the read group including PL

if [ ! -s {bam_add_rg} ];then
    picard AddOrReplaceReadGroups \\
        INPUT={bam_recal_indel} \\
        OUTPUT={bam_add_rg} \\
        RGID={tmpforrg} \\
        RGLB={tmpforrg} \\
        RGPL=ILLUMINA \\
        RGPU={tmpforrg} \\
        RGSM={tmpforrg}

    sambamba index {bam_add_rg}
else
    echo file already exist {bam_add_rg}
fi

# base recal BQSR  base quality scores recalibrate
echo -e "\\n\\ngatk BaseRecalibrator"
date

if [ ! -s {f_recal_snp_table} ];then
    # gatk3 -Xmx4g -T BaseRecalibrator -R {genome} -I {bam_recal_indel} -knownSites {dbsnp} -knownSites {mills} {known_sites_indel} -o  {f_recal_snp_table}
    gatk3 -Xmx10g -T BaseRecalibrator -R {genome} -I {bam_add_rg} -knownSites {dbsnp} -knownSites {mills} {known_sites_indel} -o {f_recal_snp_table}
else
     echo file already exist {f_recal_snp_table}
fi


echo -e "\\n\\ngatk PrintReads"
date
if [ ! -s {bam_recal_snp} ];then
    # gatk3 -Xmx2g -T PrintReads -R {genome} -I {bam_recal_indel} --BQSR {f_recal_snp_table} -o {bam_recal_snp}
    gatk3 -Xmx10g -T PrintReads -R {genome} -I {bam_add_rg} --BQSR {f_recal_snp_table} -o {bam_recal_snp}
    # gatk --java-options "Xmx10g" -T ApplyBQSR -R {genome} -I {bam_add_rg} --bqsr-recal-file {f_recal_snp_table} -o {bam_recal_snp}
else
    echo file already exist {bam_recal_snp}
fi

# MarkDuplicates
echo "picard MarkDuplicates"
date
if [ ! -s {bam_markdup} ];then
    picard -Xmx2g MarkDuplicates I={bam_recal_snp} O={bam_markdup} M={pw}/{step}/result/{lb}.markdup.metrics.txt
else
    echo file already exist {bam_markdup}
fi

# cleaning
# rm {bam_add_rg}* 2>/dev/null
# rm {bam_recal_indel}* 2>/dev/null
# rm {bam_recal_snp}* 2>/dev/null

if [ -s {bam_markdup} ];then
    echo 'completed@@@@@@@'
    touch {filedone}
fi

        """, file=out)

    # build pbs
    fl_pbs = build_pbs(pw, step, result[0], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, filedone]


def merge_bam(pw, bam_list, lb, email, dock_file):
    step = 'bam_merge'
    mds(f'{pw}/{step}', 'pbs,log,result,check')

    # get result
    result = [f'{pw}/{step}/result/{lb}.merge.bam']

    bam_list_str = '  \\\n    '.join(bam_list)

    # build pbs.sh
    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash
echo "merge sam"
date
sambamba merge {pw}/{step}/result/{lb}.merge.bam \\
    {bam_list_str}
echo "merge_complete"
date

echo "build bam index"
date
sambamba index {pw}/{step}/result/{lb}.merge.bam
date

if [ -s {result[0]} ];then
    echo 'completed@@@@@@@'
    touch {pw}/{step}/check/{lb}.{step}.done
fi
        """, file=out)

    # build pbs
    fl_pbs = build_pbs(pw, step, result[0], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done']


def call_sv(pw, bam_in, lb, ref_info, email, dock_file, tmp_root=''):
    step = 'call_sv'
    tmp_root = tmp_root if tmp_root else f'{pw}/tmp'
    pwtmp = f'{tmp_root}/call_sv'
    mds(pwtmp)
    mds(f'{pw}/{step}', 'pbs,log,result,check')

    excludebed = ref_info['excludebed']  # gm, dbsnp, mills_indel, known_indel
    genome = ref_info['gm']  # gm, dbsnp, mills_indel, known_indel
    thread = 1

    # get result
    result = [f'{pw}/{step}/result/{lb}-smoove.genotyped.vcf.gz']

    # build pbs.sh
    with open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w') as out:
        print(f"""#!/usr/bin/env bash
echo "smoove call"
date
smoove call --outdir {pw}/{step}/result --exclude {excludebed} --name {lb} --fasta {genome} -p {thread} --excludechroms "hs37d5,~:,~^GL,~decoy,~alt" --genotype {bam_in}
date
if [ -s {result[0]} ];then
    echo 'completed@@@@@@@'
    touch {pw}/{step}/check/{lb}.{step}.done
fi
        """, file=out)

    # build pbs
    fl_pbs = build_pbs(pw, step, result[0], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done']


def region_bam(step, pw, bam_merge, region_bed_file, lb, email, dock_file):
    pwstep = f'{pw}/{step}'
    mds(pwstep, 'pbs,log,result,check')

    regions = [_.strip().split('\t') for _ in open(region_bed_file) if _.strip()]
    region_refine = []
    for _ in regions:
        chr_, s, e = _[:3]

        try:
            chr_ = re.match(r'.*([\dXYM]+)', chr_).group(1)
            chr_ = {'23': 'X', '24': 'Y', '25': 'M'}[chr_] if chr_ in {'23', '24', '25'} else chr_
        except:
            print('in bed file, bad chr found:\t', _)
            continue

        try:
            region_name = _[3]
        except:
            region_name = f'{chr_}_{s}_{e}'

        region_refine.append([f'chr{chr_}:{s}-{e}', region_name])

    # print(region_refine)
    out = open(f'{pw}/{step}/pbs/{lb}_{step}.pbs.sh', 'w')
    print(f"#!/usr/bin/env bash", file=out)
    result = []
    for irg, iname in region_refine:
        ibam_result = f'{pwstep}/result/{lb}.{iname}.bam'
        print(f'samtools view -bS {bam_merge} "{irg}" -o {ibam_result} &&\\', file=out)
        print(f'sambamba index {ibam_result}', file=out)
        result.append(ibam_result)
    out.close()

    # build pbs
    fl_pbs = build_pbs(pw, step, result[-1], lb, email, dock_file, mem='40G')

    return [fl_pbs, result, f'{pw}/{step}/check/{lb}.{step}.done']



def build_pbs(pw, step, result_expected, lb, email, dock_file, mem='40G'):
    """
    for submission, we need one sh file, which would contain the actual commands to execute, this would be generated in the step function/module
    the other is a pbs file, which would examine if the result file already exist, if exist , would just ignore the run
    """
    fl_pbs = f'{pw}/{step}/pbs/{lb}_{step}.pbs'
    bash = f'singularity exec -B /mnt,/scratch {dock_file} bash ' if dock_file else 'bash'
    with open(fl_pbs, 'w') as out:
        cmd = f"""
if [ -s "{result_expected}" ];then
    echo "step {step} is already done, if need to rerun, please delete {result_expected} and then retry"
    exit 0
else
    echo "start {step} - {lb}"
    date
    {bash} {pw}/{step}/pbs/{lb}_{step}.pbs.sh &&  echo "completed"
    date
fi
"""
        pbs_cmd = pbsg(cmd, f'{pw}/{step}/log/{step}_{lb}.log', f'{step}_{lb}', email, mem=mem)
        print(pbs_cmd, file=out)

    return fl_pbs


def parse_task_dependency(task_dependency, pw, dock_file, ignore_lb=''):
    """
    parse the pbs list and generate the final submission scirpt
    each element is a task,  and consist of 2 parts,
    [task_now.pbs,  [task_depend1.pbs, task_depend2.pbs,....],  expected_result]
    ignore_lb = ignore scripts from certain steps
    """

    step = 'final_submit'
    mds(f'{pw}/{step}')

    pw_ex = f'{pw}/pbs_example'
    mds(pw_ex)
    s_step = set()

    unsubmitted = open(f'{pw}/{step}/scripts_not_submitted.txt', 'w')

    dict_sn = {}  # key=pbs name , value = sn
    n = 0
    fl_jobid = f'{pw}/{step}/jobid.list.txt'
    ignore_lb = ignore_lb if ignore_lb else []

    with open(f'{pw}/{step}/final_submit.sh', 'w') as out:
        print(f'echo start submitting >{fl_jobid}', file=out)
        for pbsnow, dep, result_exp, sample_lb in task_dependency:
            # ignore scripts from certain steps
            tmp = result_exp.rsplit('.', 2)
            if tmp[-1] != 'done':
                print('error done file', result_exp)
                continue
            step_pbs = tmp[-2]
            if step_pbs not in s_step:
                s_step.add(step_pbs)
                os.system(f'cp {pbsnow}* {pw_ex}')

            drop = 0
            for _ in ignore_lb:
                if result_exp.find(_) > -1:
                    drop = 1
                    break
            if drop:
                print(pbsnow, file=unsubmitted)
                continue

            n += 1
            dict_sn[pbsnow] = n
            #  1 -eq $1 means that , if you run this script with a argument 1, then, each step would all be submitted, no matter whether the result file exist

            dep_jid = [f'$jid{dict_sn[i]}' for i in dep if dict_sn.get(i)]
            dep_cmd = f'--dependency=afterany:{":".join(dep_jid)} ' if dep_jid else ''
            print(
                f"""

if [[ (1 -eq $1) || (! -s {result_exp}) ]]; then
    jid{n}=$(sbatch {dep_cmd} {pbsnow} | awk '{{print $NF}}')
    echo jid{n} $jid{n} {pbsnow} >> {fl_jobid}
else
    jid{n}=1000000000
fi

""", file=out)

    local_all = open(f'{pw}/{step}/final_local_run.sh', 'w')
    fh_sample = {}
    for pbsnow, dep, result_exp, sample_lb in task_dependency:
        drop = 0
        for _ in ignore_lb:
            if result_exp.find(_) > -1:
                drop = 1
                break
        if drop:
            continue


        print(f'bash {pbsnow}', file=local_all)
        try:
            print(f'bash {pbsnow}', file=fh_sample[sample_lb])
        except:
            fh_sample[sample_lb] = open(f'{pw}/{step}/local_single_sample_{sample_lb}.sh', 'w')
            print(f'bash {pbsnow}', file=fh_sample[sample_lb])
    unsubmitted.close()
    local_all.close()
    for v in fh_sample.values():
        v.close()

def generate_report():
    pass

########
# supporting functions
########


def get_read_group_from_fastq(ifl):
    """
    # @ST-E00274:188:H3JWNCCXY:4:1101:5142:1221 1:N:0:NTTGTA.
    header=$(gunzip -c $vcf|head -n 1)
    id=$(cut -f 1-4 -d: <<<$header|sed 's/@//;s/:/_/g')
    # ST-E00274_188_H3JWNCCXY_4
    sm=$(grep -Eo "[ATCGN]+$" <<<$header)
    lb=$id'_'$sm
    """
    ext = ifl.rsplit('.', 1)[-1]
    if ext == 'gz':
        header = os.popen(f'gunzip -c {ifl}|head -n 1').read().strip()
    elif ext in ['fastq', 'fq']:
        header = os.popen(f'head -n 1 {ifl}').read().strip()
    try:
        lb = re.sub('^@', '', header).split(':', 4)[:4]
        lb = '_'.join(lb)
        sm = re.match('.+([ATCGN]+$)', header).group(0)
        rg_id = f'{lb}_{sm}'
    except:
        return 'err'

    return f'@RG\\tID:{rg_id}\\tSM:{sm}\\tLB:{rg_id}\\tPU:{rg_id}\\tPL:ILLUMINA'


def getxml(srr):
    esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
        'esearch.fcgi?db=%s&term=%s&retmax=100000' % ('sra', srr)
    try:
        esearch_response = urllib.request.urlopen(esearch_url)
    except:
        return 'err'

    esearch_xml = esearch_response.read()
    esearch_root = Et.fromstring(esearch_xml)
    uids = [x.text for x in esearch_root.findall('./IdList/Id')]
    esearch_response.close()
    uids_str = ','.join(uids)

    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' + \
        f'?dbfrom=sra&db=sra&id={uids_str}'
    try:
        efetch_response = urllib.request.urlopen(efetch_url)
    except:
        time.sleep(1)
        try:
            efetch_response = urllib.request.urlopen(efetch_url)
        except:
            return 'err'

    efetch_xml = efetch_response.read()

    return efetch_xml


def get_reads_from_ncbi(srr):
    """
    get the reads/spots from ncbi
    """

    try:
        efetch_xml = os.popen(f'esearch -db sra -query "{sra_uid}"|efetch 2>/dev/null').read()
    except:
        for i in range(10):
            efetch_xml = getxml(srr)
            if efetch_xml and efetch_xml != 'err':
                break
            time.sleep(i)

    if efetch_xml and efetch_xml != 'err':
        pass
    else:
        return 'err'

    efetch_root = Et.fromstring(efetch_xml)
    try:
        run = efetch_root.find('./EXPERIMENT_PACKAGE/RUN_SET/RUN')
        return run.get('total_spots')
    except:
        return 'err'


def get_ref(genome_build):
    """
    return gm(genome_fastq), excludebed, dbsnp(vcf), mills_indel(vcf), known_indel(vcf)
    """
    ref = {}
    if genome_build == 'hg19':
        ref['gm'] = '/scratch/h_vangard_1/chenh19/ref/hg19/Homo_sapiens_assembly19.fasta'
        ref['dbsnp'] = '/scratch/h_vangard_1/chenh19/ref/gatk/dbsnp_138.b37.vcf.gz'
        ref['mills_indel'] = '/scratch/h_vangard_1/chenh19/ref/gatk/Mills_and_1000G_gold_standard.indels.b37.vcf.gz'
        ref['known_indel'] = ''
        ref['excludebed'] = '/scratch/h_vangard_1/chenh19/ref/svcall/ceph18.b37.lumpy.exclude.2014-01-15.bed'

    elif genome_build == 'hg38':
        ref['gm'] = '/scratch/h_vangard_1/chenh19/ref/hg38/Homo_sapiens_assembly38.fasta'
        ref['dbsnp'] = '/scratch/h_vangard_1/chenh19/ref/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz'
        ref['mills_indel'] = '/scratch/h_vangard_1/chenh19/ref/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
        ref['known_indel'] = '/scratch/h_vangard_1/chenh19/ref/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz'
        ref['excludebed'] = '/scratch/h_vangard_1/chenh19/ref/svcall/exclude.cnvnator_100bp.GRCh38.20170403.bed'
    return ref


def getlb(fn, rmsn=0):
    """
    if rmsn = 1, would remove the _1 / _2 (read number from the filename)
    """
    lb = re.sub(".(gz|fastq|txt|sam|bam|fastq.gz|fq.gz|tsv)$", '', fn.rsplit('/', 1)[-1])
    if rmsn:
        return re.sub(r'_(1|2)(?=[\W_])', '', lb, count=1)
    else:
        return lb


def mds(path, other_dir=None):
    '''
`f_make_dirs(path,chr_=None,other_dir=None)`
if don't want chr_ or otherdir, set them as None, or 0
1. path, the root path to make dir
2. chr_, int, generate dirs names after 1 to `int(chr_)`
e.g. if chr_=5, will make root/1, root/2, root/3.... root/5 , totally 5 folders
3. other_dir, shoud be a list,
specify the subfolders name manually
    '''
    path = f_refine_path(path)
    if not os.path.exists(path):
        os.makedirs(path)
    if other_dir:
        other_dir = newlist(other_dir)
        for i in other_dir:
            path_sub = "%s/%s" % (path, i)
            if not os.path.exists(path_sub):
                os.makedirs(path_sub)


def f_refine_path(path):
    '''
`f_refine_path(path)`
if the path is ended with backslash `/`
delete the `/`
        '''
    try:
        return re.sub('/*$', '', path)
    except:
        return ''


def newlist(l1, sep=','):
    'remove the empty elements, and strip the write space. return 0 if input datatype is not list or str'
    if isinstance(l1, str):
        l1 = l1.split(sep)
    elif not isinstance(l1, list):
        return 0
    return [i.strip() if isinstance(i, str) else i for i in l1 if i]


def pbsg(cmd, log, lb, email, mem="40G", timeout=72, nnode=1):
    """return the format for pbs script, 1 = cmd, 2 = log, 3=email, 4=max-mem, 5 = timeout, default 72h"""

    job = f'#SBATCH --job-name={lb}' if lb else ''
    t = timeout if timeout not in [0, 'na', "NA"] else 72
    info = """#!/bin/bash
#SBATCH --mail-user=%s
#SBATCH --mail-type=ALL
#SBATCH --nodes=%s
#SBATCH --ntasks=1
#SBATCH --time=%s:00:00
#SBATCH --mem=%s
#SBATCH -o %s
%s

%s

""" % (email, nnode, t, mem, log, job, cmd)
    return info
