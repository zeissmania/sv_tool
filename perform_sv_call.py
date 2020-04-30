"""
parse the sv call config file
"""
import sys, os
# load logging
import logging
import logging.config
import yaml
from time import time as t
sys.path.append(f'/pipeline')   # withiin singularity image
from sv_call import sv_call_modules

# steps: ['sra2fastq ', 'split_fastq', 'bwa_align', 'bam_refine', 'bam_merge', 'call_sv']
# if you don't want to submit scripts from certain step, just add them here

toggle_time = 0
time_count = {_:0 for _ in ['init', 'pos1', 'pos2', 'pos3', 'ftp', 'sra2fastq', 'split_fastq', 'read_group', 'bwa_align', 'bam_refine', 'bam_merge', 'call_sv', 'parse_pbs']}

def run(config, config_default):

    s = t()
    # read the arguments
    pw = config.get('pw')
    assert pw is not None, 'pw, working path not set ind config'
    email = config.get('email')
    assert email is not None, 'email not found in config'
    gm_build = config.get('gm_build')
    assert gm_build is not None, 'gm_build not found in config'
    sra_convert = config.get('sra_to_fastq')
    assert sra_convert is not None, 'sra_to_fastq not found in config'
    download_sra = config.get('download_sra')
    assert download_sra is not None, 'download_sra not found in config'
    filelist = config.get('filelist')
    assert filelist is not None, 'filelist not found in config'
    paired = config.get('paired')
    assert paired is not None, 'paired not found in config'

    read_group_fq = config_default.get('read_group_from_fastq')
    assert read_group_fq is not None, 'read_group_from_fastq not found in config_default'
    chunk = config_default.get('chunk')
    assert chunk is not None, 'chunk not found in config_default'
    dock_file = config_default.get('dock_file')
    assert dock_file is not None, 'dock_file not found in config_default'
    no_cal = config_default.get('no_cal')
    assert no_cal is not None, 'no_cal not found in config_default'
    tmp_dir = config_default.get('tmp_dir')
    assert tmp_dir is not None, 'tmp_dir is not found in config_default'

    ignore_lb = config_default.get('ignore_lb')
    assert ignore_lb is not None, 'ignore_lb not found in config_default'

    region_bed_file = config_default.get('region_bed_file')

    # if compress_when_split = True, would compress the file when spliting,
    # this would significantly reduce the max tmp file size
    # but would be slower if running under accre cluster
    # because when under accre cluster  the compression is performed at the bwa_map step, would run in parallel tasks
    # the gzip step  is sequential when this toggle is true
    compress_when_split = True


    fl_log_conf = f'/home/chenh19/jb/config/logging_setting.yaml'
    log_prefix = 'build_sv_pbs'
    log_prefix = log_prefix + '_' if log_prefix else ''
    with open(fl_log_conf) as f:
        cfg = yaml.safe_load(f.read())
        cfg['handlers']['file']['filename'] = log_prefix + cfg['handlers']['file']['filename']
        cfg['handlers']['error']['filename'] = log_prefix + cfg['handlers']['error']['filename']
        logging.config.dictConfig(cfg)
    logger = logging.getLogger('main')

    logger.debug('\n'*3+'*'*20+'\n\tnew_run\n' + "*"*20)

    if tmp_dir:
        if os.path.isdir(tmp_dir):
            pass
        else:
            try:
                os.mkdir(tmp_dir)
            except:
                logger.warning(f'the tmpdir in config cannot be build, the tmp dir would be written to  {pw}/tmp')
                tmp_dir = f'{pw}/tmp'

    logger.debug(f'the split tmp file dir = {tmp_dir}')


    # if we need to convert sra, we can't get readgroup from fastq
    if sra_convert:
        read_group_fq = 0

    ref_info = sv_call_modules.get_ref(gm_build)  # return = dict . gm(genome_fastq), exclude_bed, dbsnp(vcf), mills_indel(vcf), known_indel(vcf)
    genome = ref_info['gm']

    err = sv_call_modules.validate_sv_config(pw, gm_build, sra_convert, filelist, download_sra, dock_file, ref_info, config_default)
    if err:
        print('pbs is not generated due to the following error')
        print('\n'.join(err))
        return 1

    task_dependency = []
    # each element is a task,  and consist of [task_now.pbs,  [task_depend1.pbs, task_depend2.pbs,....]]


    if toggle_time:
        step = 'init'
        time_count[step] += t() -s
        s = t()

    # if some sample list has error, would ignore these scripts when submitting

    for lb, ifl in filelist.items():

        # replace all the dot in lb, because in the following steps, we may need dot as sep to extract info
        lb = lb.replace('.', '_')

        # sra_convert
        if sra_convert:
            # sra2fastq_dryrun = 1 if 'sra2fastq' in ignore_lb else 0
            try:
                tmp = sv_call_modules.sra2fastq(pw, lb, ifl[0], paired, email, dock_file)
                pbs_sra_convert, fastq_files, donefile, time_count1 = tmp

            except:
                # the return of sra2fastq is error. skip the following steps
                continue
            time_count['pos1'] += time_count1['pos1']
            time_count['pos2'] += time_count1['pos2']
            time_count['pos3'] += time_count1['pos3']
            time_count['ftp'] += time_count1['ftp']
        # fastq_files = [f'{pw}/sra2fastq/result/{sra_file}_1.fastq.gz',
        #       f'{pw}/sra2fastq/result/{sra_file}_2.fastq.gz'] if paired else [f'{pw}/sra2fastq/result/{sra_file}.fastq.gz']
            task_dependency.append([pbs_sra_convert, [], donefile, lb])
        else:
            fastq_files = ifl
            pbs_sra_convert = ''

        if toggle_time:
            step = 'sra2fastq'
            time_count[step] += t() -s
            s = t()

        # split fastq
        pbs_split_fastq, fastq_files_set, donefile = sv_call_modules.split_fastq(pw, lb, fastq_files, email, dock_file, chunk=chunk, tmp_root=tmp_dir, compress=compress_when_split)
        # f'{pw}/tmp'
        task_dependency.append([pbs_split_fastq, [pbs_sra_convert] if pbs_sra_convert else [], donefile, lb])
        # result = [(task_1_chunkx001.fastq.gz, task_2_chunkx001.fastq.gz), (task_1_chunkx002.fastq.gz, task_2_chunkx002.fastq.gz), ...]

        if toggle_time:
            step = 'split_fastq'
            time_count[step] += t() -s
            s = t()


        # get read group
        read_group_raw = f'@RG\\tID:{lb}\\tSM:{lb}\\tLB:{lb}\\tPU:{lb}\\tPL:ILLUMINA'

        if read_group_fq:
            read_group = sv_call_modules.get_read_group_from_fastq(fastq_files[0])
            if read_group == 'err':
                # print(f'fail to get read_group from fastq, use key from filelist as readgroup: {lb} : {ifl}')
                read_group = read_group_raw
        if not read_group_fq:
            read_group = read_group_raw

        if toggle_time:
            step = 'read_group'
            time_count[step] += t() -s
            s = t()

        bam_chunks = []
        pbs_recal_all = []
        for fastq_chunk in fastq_files_set:
            # bwa
            pbs_bwa, bam_bwa, donefile = sv_call_modules.bwa_align(pw, fastq_chunk, genome, email, read_group, dock_file)
            # bam_bwa = [f'{pw}/{step}/result/{lb}.sorted.bam']
            bam_bwa = bam_bwa[0]
            task_dependency.append([pbs_bwa, [pbs_split_fastq], donefile, lb])

            if toggle_time:
                step = 'bwa_align'
                time_count[step] += t() -s
                s = t()

            # recal
            if no_cal:
                pbs_recal_all.append(pbs_bwa)
                bam_chunks.append(bam_bwa)
            else:
                pbs_recal, bam_recal, donefile = sv_call_modules.refine_bam(pw, bam_bwa, ref_info, email, dock_file, lb)
                # bam_recal = [f'{pw}/{step}/result/lb.recal.markdup.bam']
                bam_recal = bam_recal[0]
                task_dependency.append([pbs_recal, [pbs_bwa], donefile, lb])

                pbs_recal_all.append(pbs_recal)
                bam_chunks.append(bam_recal)

            if toggle_time:
                step = 'bam_refine'
                time_count[step] += t() -s
                s = t()

        # merge
        pbs_merge, bam_merge, donefile = sv_call_modules.merge_bam(pw, bam_chunks, lb, email, dock_file)
        # bam_merge = [f'{pw}/{step}/result/{lb}.merge.bam']
        bam_merge = bam_merge[0]
        task_dependency.append([pbs_merge, pbs_recal_all, donefile, lb])

        if toggle_time:
            step = 'bam_merge'
            time_count[step] += t() -s
            s = t()

        # call sv
        pbs_callsv, vcf_call_sv, donefile = sv_call_modules.call_sv(pw, bam_merge, lb, ref_info, email, dock_file)
        # vcf_call_sv = [f'{pw}/{step}/result/{lb}-smoove.genotyped.vcf.gz']
        vcf_call_sv = vcf_call_sv[0]
        task_dependency.append([pbs_callsv, [pbs_merge], donefile, lb])

        if toggle_time:
            step = 'call_sv'
            time_count[step] += t() -s
            s = t()

        # extract the regional bam file

        if region_bed_file:
            step = 'extract_region_bam'
            pbs_region_bam, bam_region, donefile = sv_call_modules.region_bam(step, pw, bam_merge, region_bed_file, lb, email, dock_file)
            task_dependency.append([pbs_region_bam, [pbs_merge], donefile, lb])


    # parse the task_dependency, generate the submission task
    sv_call_modules.parse_task_dependency(task_dependency, pw, dock_file, ignore_lb)

    if toggle_time:
        step = 'parse_pbs'
        time_count[step] += t() -s
        s = t()
        for k, v in time_count.items():
            print(f'{k}\t{round(v,2)}')
