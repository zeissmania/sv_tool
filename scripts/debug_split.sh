#single end
python ../debug.py split_fastq -i ~/test/download/SRR4253621.fastq.gz -o ~/test/split/SRR4253621 --trunk 10

#pair end
python ../debug.py split_fastq -i ~/test/download/SRR4032155_1.fastq.gz,~/test/download/SRR4032155_2.fastq.gz -o ~/test/split/SRR4032155 --trunk 10
