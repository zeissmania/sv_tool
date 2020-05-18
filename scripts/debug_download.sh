if [[! -s /scratch/cqs/temp/download]]; then
  mkdir /scratch/cqs/temp/download
fi

#single end
python ../debug.py download -i SRR4253621 -s -o /scratch/cqs/temp/download

#pair end
python ../debug.py download -i SRR4032155 -s -o /scratch/cqs/temp/download
