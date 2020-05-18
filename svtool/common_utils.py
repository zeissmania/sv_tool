import os
import os.path
import logging
import errno
import shutil

from .__version__ import __version__, __name__

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def run_command(logger, cmd):
  logger.info(cmd)
  os.system(cmd)

def read_file_map(fileName):
  check_file_exists(fileName)

  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def check_file_map(fileMap):
  for sname in fileMap.keys():
    sfile = fileMap[sname]
    check_file_exists(sfile)

def initialize_logger(logfile, name, args, loglevel=logging.DEBUG):
  logger = logging.getLogger(name)
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  logger.info(f"{__name__}: version={__version__}, command={name}, args={args}")
  return(logger)

def write_r_script(outfilePrefix, rScript, optionMap={}):
  check_file_exists(rScript)

  targetScript = outfilePrefix + ".r"
  optionMap["outfilePrefix"] = outfilePrefix
  with open(targetScript, "wt") as fout:
    for key in optionMap.keys():
      fout.write("%s='%s'\n" % (key, optionMap[key]))

    fout.write("setwd('%s')\n" % os.path.dirname(os.path.abspath(targetScript)))
    
    with open(rScript, "rt") as fin:
      bFirstSetwd = True
      for line in fin:
        if line.startswith("setwd") and bFirstSetwd:
          bFirstSetwd = False
          continue

        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "="):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)

def write_rmd_script(outfilePrefix, rmdScript, optionMap={}, copyRFunction=True, optionsToIndividualFile=True):
  check_file_exists(rmdScript)

  if copyRFunction:
    rFunScript = os.path.join( os.path.dirname(__file__), "Rfunctions.R")
    check_file_exists(rFunScript)

    targetFolder = os.path.dirname(os.path.abspath(outfilePrefix))
    targetRFunScript =  os.path.join(targetFolder, "Rfunctions.R")
    shutil.copyfile(rFunScript, targetRFunScript)

  if optionsToIndividualFile:
    optionToIndividualFileName = outfilePrefix + ".options"
    with open(optionToIndividualFileName, "wt") as fout:
      for key in sorted(optionMap.keys()):
        fout.write("%s\t%s\n" % (key, optionMap[key]))
    optionMap = {"option_file": os.path.basename(optionToIndividualFileName)}

  targetScript = outfilePrefix + ".rmd"
  with open(targetScript, "wt") as fout:
    with open(rmdScript, "rt") as fin:
      for line in fin:
        if line.startswith("```"):
          fout.write(line)
          for key in optionMap.keys():
            fout.write("%s='%s'\n" % (key, optionMap[key]))
          fout.write("\n")
          break
        else:
          fout.write(line)

      for line in fin:
        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "=") or line.startswith(key + "<-"):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)
