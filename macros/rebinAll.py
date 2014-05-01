#!/usr/bin/env python
"""Rebin all histograms in a file"""
__version__ = "1.0"

import sys
import optparse
import shutil
import os
import re

if '-h' not in sys.argv and len(sys.argv) > 0:
    import ROOT
    # ROOT parses options when the first ROOT command is called, so we must
    # add '-b' before that to get batch mode, but we must immediately remove
    # it to avoid interference with option parsing for this script.
    sys.argv.append('-b')
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    sys.argv.remove('-b')

class RootFile:
    """A wrapper for TFiles, allowing quick access to the name and Get."""
    def __init__(self, file_name):
        self.name = file_name[0:-5]
        self.file = ROOT.TFile(file_name, "read")
        if self.file.IsZombie():
            print "Error opening %s, exiting..." % file_name
            sys.exit(1)
    def Get(self, object_name):
        return self.file.Get(object_name)

def process_directory(path,files):
    """Loop through all histograms in the directory and plot them."""
    newfile = ROOT.TFile(options.newfile, "RECREATE")
    keys = files[0].file.GetDirectory(path).GetListOfKeys()
    key = keys[0]
    while key:
        obj = key.ReadObj()
        key = keys.After(key)
        if (obj.IsA().InheritsFrom("TH1") and
            not obj.IsA().InheritsFrom("TH2") and
            not obj.IsA().InheritsFrom("TH3")):
            newfile.cd()
            newobj = obj.Clone()
            newobj.Rebin(options.nbins)
            newobj.GetXaxis().SetRangeUser(900,2500)
    newfile.Write()
    newfile.Close()

def main():
  usage="""Usage: %prog [options] file.root"""
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-b', '--bins', default="1", metavar="NUMBER", type="int", help="Number of bins to merge")
  parser.add_option('--output', default="rebinned.root", metavar="NAME", help="Name of rebinned file; default is 'rebinned.root'")
  global options
  options, arguments = parser.parse_args()
  options.newfile = options.output
  options.nbins = options.bins
  print options.newfile, options.nbins
  files = [RootFile(filename) for filename in arguments]
  process_directory("/", files)


if __name__ == "__main__":
  main() 
