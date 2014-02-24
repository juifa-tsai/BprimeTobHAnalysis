#!/usr/bin/python
import os, sys, math, ROOT

ROOT.gROOT.SetBatch() 

if __name__  == "__main__":
  try: 
    rootFileName = sys.argv[1]
  except IndexError:
    sys.stderr.write("Missing root file argument.")
    sys.stderr.write(os.linesep)
    sys.exit(-1)

  if not os.path.exists(rootFileName):
    sys.stderr.write("Incorrect root file argument \"%s\"." % rootFileName)
    sys.stderr.write(os.linesep)
    sys.exit(-1)

  file = ROOT.TFile(rootFileName)

  hdata = file.Get("DATA__h_cutflow")
  htt   = file.Get("TTJets__h_cutflow") ; 
  hqcd  = file.Get("QCD__h_cutflow") ; 
  hsig = [ 
   file.Get("BprimeToBHinc_M-500_8TeV-madgraph__h_cutflow") , 
   file.Get("BprimeToBHinc_M-600_8TeV-madgraph__h_cutflow") , 
   file.Get("BprimeToBHinc_M-800_8TeV-madgraph__h_cutflow") , 
   file.Get("BprimeToBHinc_M-1000_8TeV-madgraph__h_cutflow"),  
   file.Get("BprimeBprimeToBHBHinc_M-500__h_cutflow") , 
   file.Get("BprimeBprimeToBHBHinc_M-600__h_cutflow") , 
   file.Get("BprimeBprimeToBHBHinc_M-700__h_cutflow") , 
   file.Get("BprimeBprimeToBHBHinc_M-800__h_cutflow") , 
   file.Get("BprimeBprimeToBHBHinc_M-1000__h_cutflow"),  
   file.Get("BprimeBprimeToBHBHinc_M-1200__h_cutflow")   
  ] 

  hmc = htt.Clone("hmc")
  hmc.Add(hqcd)

  for hist in hsig:
    fname = 'datacard_' + hist.GetName().split('__')[0] + '.txt'
    mass =  (hist.GetName().split('M-')[1]).split('_')[0]
    outf = open(fname, 'w')
    outf.write("-------------------------------------------------------------------------------------------------------------------\n")
    outf.write("imax 2\n")    
    outf.write("jmax 2\n")    
    outf.write("kmax *\n")    
    outf.write("------------------------------------------------------------------------------------\n")  
    outf.write("bin         1b     2b\n") 
    outf.write("observation " + str(hdata.GetBinContent(8)) + " " + str(hdata.GetBinContent(9)) + "\n") 
    outf.write("------------------------------------------------------------------------------------\n")  
    outf.write("bin                1b         1b         1b         2b         2b         2b\n") 
    outf.write("process            bp"+str(mass)+"      ttjets     qcd        bp"+str(mass)+"      ttjets     qcd\n")  
    outf.write("process            0          1          2          0          1          2\n")  
    outf.write("rate               " + str(round(hist.GetBinContent(8), 3))+ "     "+ str(round(htt.GetBinContent(8), 3))+ "     "+ str(round(hqcd.GetBinContent(8), 3))+  "    " + str(round(hist.GetBinContent(9), 3))+ "       "+ str(round(htt.GetBinContent(9), 3))+ "     "+ str(round(hqcd.GetBinContent(9), 3))+"\n")
    outf.write("------------------------------------------------------------------------------------\n")  
    outf.write("lumi        lnN    1.026      1.026      1.026      1.026      1.026      1.026\n")
    outf.write("norm_sig    lnN    1.30       -          -          1.30       -          -\n") 
    outf.write("norm_bkg    lnN    -          1.30       1.30       -          1.30       1.30\n") 
    outf.write("-------------------------------------------------------------------------------------------------------------------\n")
    outf.close() 

