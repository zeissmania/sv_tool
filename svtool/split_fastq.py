import argparse
import sys
import logging
import os
import gzip
import math
from _ctypes import ArgumentError

from .common_utils import check_file_exists, initialize_logger

def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)

def rawgencount(filename):
    with gzip.open(filename, 'rb') as f:
      f_gen = _make_gen(f.read)
      return sum( buf.count(b'\n') for buf in f_gen )

def split_fastq(logger, inputFileStr, outputFilePrefix, trunkNumber):  
  inputFiles = inputFileStr.split(",")
  
  if len(inputFiles) == 1:
    isPairedEnd = False
  elif len(inputFiles) == 2:
    isPairedEnd = True
  else:
    raise ArgumentError('inputFile should be only one file (single end) or two files (pair end): %s ' % inputFileStr)

  for inputFile in inputFiles:
    check_file_exists(inputFile)
    
  fileIndex = 0
  for inputFile in inputFiles:
    fileIndex = fileIndex + 1
    logger.info("Processing %s ..." % inputFile)
    totalLineCount = rawgencount(inputFile)
    totalRecord = totalLineCount / 4
    recordPerFile = math.ceil(totalRecord / trunkNumber)
    logger.info("Total %d reads, each file should have almost %d reads." % (totalRecord, recordPerFile))
    
    with gzip.open(inputFile, "rt") as fin:
      for trunk in range(1, (trunkNumber + 1)):
        if isPairedEnd:
          fileName = "%s.%d.%d.fastq.gz" % (outputFilePrefix, trunk, fileIndex)
        else:
          fileName = "%s.%d.fastq.gz" % (outputFilePrefix, trunk)
        
        logger.info("Writing reads to %s ..." % fileName)
        with gzip.open(fileName, "wt") as fout:
          for idx in range(0, recordPerFile):
            line1 = fin.readline()
            if not line1:
              break
            fout.write(line1)
            fout.write(fin.readline())
            fout.write(fin.readline())
            fout.write(fin.readline())
        
  logger.info("done")

def add_split_parser(parser, NOT_DEBUG=True):
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)
  parser.add_argument('--trunk', action='store', nargs='?', type=int, default=50, help="Number of small files")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Split big fastq file to multiple small fastq files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  add_split_parser(parser, NOT_DEBUG)

  args = parser.parse_args()
  
  if DEBUG:
    args.input = "W:/SequenceData/20191121_4145-DM-1/4145-DM-1_1-AACCAGAT-TCTTTCCC_S99_R1_001.fastq.gz,W:/SequenceData/20191121_4145-DM-1/4145-DM-1_1-AACCAGAT-TCTTTCCC_S99_R2_001.fastq.gz"
    args.outputPrefix = "E:/temp/splitFastqTest"
    args.trunk = 3
  
  logger = logging.getLogger('splitFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  split_fastq(logger, args.input, args.outputPrefix, args.trunk)
  
if __name__ == "__main__":
    main()
