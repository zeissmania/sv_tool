#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $target_folder = "/scratch/cqs/temp";

my $config = {
  general => {
    task_name => "svtool",
  },
  "download_single_end" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_folder}/download_single_end",
    option                => "download -s", #single end
    interpretor           => "",
    check_program         => 0,
    program               => "svtool",
    source                => {"SRR4253621" => ["SRR4253621"]},
    source_arg            => "-i",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_ext       => ".fastq.gz",
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "bwa_single_end" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_folder}/bwa_single_end",
    option                => "",
    interpretor           => "",
    check_program         => 0,
    program               => "bwa",
    source_ref            => "download_single_end",
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".bam",
    output_file_ext       => ".bam",
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "download_pair_end" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_folder}/download_pair_end",
    option                => "download", #pair end
    interpretor           => "",
    check_program         => 0,
    program               => "svtool",
    source                => {"SRR4032155" => ["SRR4032155"]},
    source_arg            => "-i",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_ext       => "_1.fastq.gz;_2.fastq.gz",
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "bwa_pair_end" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_folder}/bwa_pair_end",
    option                => "",
    interpretor           => "",
    check_program         => 0,
    program               => "bwa",
    source_ref            => "download_pair_end",
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".bam",
    output_file_ext       => ".bam",
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTaskSlurmSlim",
    perform    => 1,
    target_dir => "${target_folder}/sequencetask",
    option     => "",
    source     => {
      step1 => ["download_single_end", "download_pair_end", "bwa_single_end", "bwa_pair_end"],
    },
    sh_direct => 0,
    cluster   => "slurm",
    pbs       => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  }
};

#performTask($config, "download_single_end");
performConfig($config);

1;
