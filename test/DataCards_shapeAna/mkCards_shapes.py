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

  hdata1b = file.Get("DATA_HT_1bjet_Observe")
  hdata2b = file.Get("DATA_HT_over1bjet_Observe")
  hbkg1b = file.Get("DATA_HT_1bjet_Expect")
  hbkg2b = file.Get("DATA_HT_over1bjet_Expect")
  hsig1b = [ 
   file.Get("BHBH500_HT_1bjet_Signal"), 
   file.Get("BHBH600_HT_1bjet_Signal"), 
   file.Get("BHBH700_HT_1bjet_Signal"), 
   file.Get("BHBH800_HT_1bjet_Signal"), 
   file.Get("BHBH1000_HT_1bjet_Signal"), 
   file.Get("BHBH1200_HT_1bjet_Signal"), 
  ] 
  hsig2b = [ 
   file.Get("BHBH500_HT_over1bjet_Signal"), 
   file.Get("BHBH600_HT_over1bjet_Signal"), 
   file.Get("BHBH700_HT_over1bjet_Signal"), 
   file.Get("BHBH800_HT_over1bjet_Signal"), 
   file.Get("BHBH1000_HT_over1bjet_Signal"), 
   file.Get("BHBH1200_HT_over1bjet_Signal"), 
  ] 

  for hist1b in hsig1b:
    index = hsig1b.index(hist1b) 
    hist2b = hsig2b[index]
    mass =  (hist1b.GetName().split('BHBH')[1]).split('_')[0]
    fname = 'datacard_BHBH'+'_M-'+str(mass)+'_normalizedShapes.txt'
    print index, fname
    outf = open(fname, 'w')
    outf.write("------------------------------------------------------------------------------------\n")
    outf.write("imax 2 number of channels\n")
    outf.write("jmax 2 number of backgrounds\n")
    outf.write("kmax * number of nuisance parameters (sources of systematical uncertainties)\n")
    outf.write("------------------------------------------------------------------------------------\n")
    outf.write("shapes BHBH"+str(mass)+"  Cat1b  "+rootFileName+"  BHBH"+str(mass)+"_HT_1bjet_Signal\n")
    outf.write("shapes Bkg      Cat1b  "+rootFileName+"  DATA_HT_1bjet_Expect\n")
    outf.write("shapes data_obs Cat1b  "+rootFileName+"  DATA_HT_1bjet_Observe\n")
    outf.write("shapes BHBH"+str(mass)+"  Cat2b  "+rootFileName+"  BHBH"+str(mass)+"_HT_over1bjet_Signal\n")
    outf.write("shapes Bkg      Cat2b  "+rootFileName+"  DATA_HT_over1bjet_Expect\n")
    outf.write("shapes data_obs Cat2b  "+rootFileName+"  DATA_HT_over1bjet_Observe\n")
    outf.write("------------------------------------------------------------------------------------\n")
    outf.write("bin         Cat1b  Cat2b\n")
    outf.write("observation "+str(hdata1b.Integral())+"   "+str(hdata2b.Integral())+"\n")
    outf.write("------------------------------------------------------------------------------------\n")
    outf.write("bin                Cat1b      Cat1b      Cat1b      Cat2b      Cat2b  Cat2b\n")
    outf.write("process            BHBH"+str(mass)+"    Bkg        data_obs   BHBH"+str(mass)+"    Bkg    data_obs\n")
    outf.write("process            0          1          2          0          1      2\n")
    outf.write("rate               "+str(round(hist1b.Integral(),3))+"     "+str(round(hbkg1b.Integral(),2))+"    "+str(round(hdata1b.Integral(),3))+"     "+str(round(hist2b.Integral(),3))+"     "+str(round(hbkg2b.Integral(),3))+"  "+str(round(hdata2b.Integral(),3))+" \n") 
    outf.write("------------------------------------------------------------------------------------\n")
    outf.write("lumi        lnN    1.026      1.026      -          1.026      1.026  -\n")
    outf.write("norm_sig    lnN    1.30       -          -          1.30       -      -\n")
    outf.write("norm_bkg    lnN    -          1.30       -          -          1.30   -\n")
    outf.write("------------------------------------------------------------------------------------\n")
    outf.close() 

