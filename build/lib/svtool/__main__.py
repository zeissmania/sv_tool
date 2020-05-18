import argparse
import sys
import os
from pathlib import Path

from .__version__ import __version__
from .common_utils import initialize_logger
from .download_fastq import add_download_parser, download

def main():
  parser = argparse.ArgumentParser(description="Structure Variant Tool " + __version__,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  subparsers = parser.add_subparsers(dest="command")

  parser_d = subparsers.add_parser('download', help='Download FASTQ files based on SRR ID')
  add_download_parser(parser_d, NOT_DEBUG)
  
  if not DEBUG and len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  
  if args.command == "download":
    logger = initialize_logger(os.path.join(args.output, args.input + ".log"), args.command, args)
    download(logger, args.input, args.is_single_end, args.output)
  
if __name__ == "__main__":
    main()
