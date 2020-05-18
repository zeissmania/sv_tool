if [[! -s /scratch/cqs/temp/split]]; then
  mkdir /scratch/cqs/temp/split
fi

#single end
python ../debug.py split_fastq -i /scratch/cqs/temp/download/SRR4253621.fastq.gz -o /scratch/cqs/temp/split/SRR4253621 --trunk 10

#pair end
python ../debug.py split_fastq -i /scratch/cqs/temp/download/SRR4032155_1.fastq.gz,/scratch/cqs/temp/download/SRR4032155_2.fastq.gz -o /scratch/cqs/temp/split/SRR4032155 --trunk 10
