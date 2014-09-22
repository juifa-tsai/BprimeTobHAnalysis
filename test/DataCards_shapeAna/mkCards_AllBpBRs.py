#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT

datacard_template = """--------------------------------------------------------------------------------------------------------------------

imax 2 number of channels
jmax 1 number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
--------------------------------------------------------------------------------------------------------------------
shapes * * ROOTFILENAME $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC
--------------------------------------------------------------------------------------------------------------------
bin                       Cat1b                         Cat2b
observation               NOBSCAT1B                     NOBSCAT2B
--------------------------------------------------------------------------------------------------------------------
bin                       Cat1b                         Cat1b         Cat2b                          Cat2b     
process                   BpMMBP_bHBH_bZBZ_tWTW        background    BpMMBP_bHBH_bZBZ_tWTW         background
process                   0                             1             0                              1         
rate                      NSIGCAT1B                         NBKGCAT1B       NSIGCAT2B                          NBKGCAT2B
--------------------------------------------------------------------------------------------------------------------
lumi               lnN    1.026                         -             1.026                          -         
purewt             lnN    0.990/1.010                   -             0.990/1.010                    -         
pdfrewt            lnN    0.990/1.010                   -             0.990/1.010                    -         
trigsf             lnN    1.01                          -             1.02                           -         
ca8misc            lnN    1.15                          -             1.15                           -
JES                shape  1                             1             1                              1 
JER                shape  1                             1             1                              1 
CA8                shape  1                             1             1                              1 
SFb                shape  1                             1             1                              1 
SFl                shape  1                             1             1                              1 
SigStat            shape  1                             -             1                              -
Stat               shape  -                             1             -                              1
TTJetsScale        shape  -                             1             -                              1
TTJetsMatching     shape  -                             1             -                              1
TopPtReWrt         lnN    -                             0.992/1.005   -                              0.996/1.002
--------------------------------------------------------------------------------------------------------------------
"""

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

  hdata1b = file.Get("Cat1b_data_obs")
  hdata2b = file.Get("Cat2b_data_obs")
  hbkg1b = file.Get("Cat1b_background")
  hbkg2b = file.Get("Cat2b_background")

  counts = 0
  for BRbh in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    BRbz = 0
    while BRbz <= (10 - BRbh)+0.1:
      BRtw = 10 - BRbz - BRbh
      print ' BRbh ', BRbh, ' BRbz ', BRbz, ' BRtw ', BRtw 
      hsig1b = [ 
        file.Get("Cat1b_BpM500_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM600_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM700_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM800_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM900_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM1000_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat1b_BpM1200_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
       ] 
      hsig2b = [ 
        file.Get("Cat2b_BpM500_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM600_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM700_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM800_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM900_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM1000_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        file.Get("Cat2b_BpM1200_bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)), 
        ] 

      for hist1b in hsig1b:
        index = hsig1b.index(hist1b) 
        hist2b = hsig2b[index]
        mass =  (hist1b.GetName().split('_')[1]).split('BpM')[1]
        card = open(os.path.join('datacard_BprimebH_BpM'+str(mass)+'_bH'+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2)+'.txt'), 'w')
        card_content = re.sub('ROOTFILENAME',rootFileName,datacard_template)
        card_content = re.sub('BpMMBP','BpM'+str(mass),card_content)
        card_content = re.sub('bHBH_bZBZ_tWTW',"bH"+str(BRbh).zfill(2)+"_bZ"+str(BRbz).zfill(2)+"_tW"+str(BRtw).zfill(2),card_content)
        card_content = re.sub('NOBSCAT1B',str(round(hdata1b.Integral(),3)),card_content)
        card_content = re.sub('NOBSCAT2B',str(round(hdata2b.Integral(),3)),card_content)
        card_content = re.sub('NSIGCAT1B',str(round(hist1b.Integral(),3)),card_content)
        card_content = re.sub('NSIGCAT2B',str(round(hist2b.Integral(),3)),card_content)
        card_content = re.sub('NBKGCAT1B',str(round(hbkg1b.Integral(),3)),card_content)
        card_content = re.sub('NBKGCAT2B',str(round(hbkg2b.Integral(),3)),card_content)
        card.write(card_content)
        card.close() 

      BRbz += 1 
      counts += 1
  print counts 

