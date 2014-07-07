#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT

datacard_template_Cat1b = """---------------------------------------------------------------------------------------------
imax 1 number of channels
jmax 1 number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
---------------------------------------------------------------------------------------------
shapes * * ROOTFILENAME $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------
bin                Cat1b           
observation        NOBSCAT1B
---------------------------------------------
bin                Cat1b           Cat1b     
process            BHBHMBP         background
process            0               1         
rate               NSIGCAT1B       NBKGCAT1B 
---------------------------------------------
lumi        lnN    1.026           -         
purewt      lnN    0.990/1.010     -         
pdfrewt     lnN    0.990/1.010     -         
trigsf      lnN    1.01            -         
Syst        shape  -               1         
JES         shape  1               -         
JER         shape  1               -         
SFCA8       shape  1               -         
SFHb        shape  1               -         
SFHl        shape  1               -         
SFb         shape  1               -         
SFl         shape  1               -         
---------------------------------------------
"""

datacard_template_Cat2b = """---------------------------------------------------------------------------------------------
imax 1 number of channels
jmax 1 number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
---------------------------------------------------------------------------------------------
shapes * * ROOTFILENAME $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------
bin                Cat2b
observation        NOBSCAT2B
---------------------------------------------------------------
bin                Cat2b           Cat2b     
process            BHBHMBP         background
process            0               1         
rate               NSIGCAT2B       NBKGCAT2B
----------------------------------------------------------------
lumi        lnN    1.026           -         
purewt      lnN    0.990/1.010     -         
pdfrewt     lnN    0.990/1.010     -         
trigsf      lnN    1.02            -         
Syst        shape  -               1
JES         shape  1               - 
JER         shape  1               - 
SFCA8       shape  1               -         
SFHb        shape  1               - 
SFHl        shape  1               - 
SFb         shape  1               - 
SFl         shape  1               - 
---------------------------------------------------------------
"""

datacard_template = """---------------------------------------------------------------------------------------------
imax 2 number of channels
jmax 1 number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
---------------------------------------------------------------------------------------------
shapes * * ROOTFILENAME $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------
bin                Cat1b                         Cat2b
observation        NOBSCAT1B                     NOBSCAT2B
---------------------------------------------------------------------------------------------
bin                Cat1b           Cat1b         Cat2b           Cat2b     
process            BHBHMBP         background    BHBHMBP         background
process            0               1             0               1         
rate               NSIGCAT1B       NBKGCAT1B     NSIGCAT2B       NBKGCAT2B
----------------------------------------------------------------------------------------------
lumi        lnN    1.026           -             1.026           -         
purewt      lnN    0.990/1.010     -             0.990/1.010     -         
pdfrewt     lnN    0.990/1.010     -             0.990/1.010     -         
trigsf      lnN    1.01            -             1.02            -         
Syst        shape  -               1             -               1
JES         shape  1               -             1               - 
JER         shape  1               -             1               - 
SFCA8       shape  1               -             1               - 
SFHb        shape  1               -             1               - 
SFHl        shape  1               -             1               - 
SFb         shape  1               -             1               - 
SFl         shape  1               -             1               - 
---------------------------------------------------------------------------------------------
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
  hsig1b = [ 
   file.Get("Cat1b_BHBH500"), 
   file.Get("Cat1b_BHBH600"), 
   file.Get("Cat1b_BHBH700"), 
   file.Get("Cat1b_BHBH800"), 
   file.Get("Cat1b_BHBH900"), 
   file.Get("Cat1b_BHBH1000"), 
   file.Get("Cat1b_BHBH1200"), 
  ] 
  hsig2b = [ 
   file.Get("Cat2b_BHBH500"), 
   file.Get("Cat2b_BHBH600"), 
   file.Get("Cat2b_BHBH700"), 
   file.Get("Cat2b_BHBH800"), 
   file.Get("Cat2b_BHBH900"), 
   file.Get("Cat2b_BHBH1000"), 
   file.Get("Cat2b_BHBH1200"), 
  ] 

  for hist1b in hsig1b:
    index = hsig1b.index(hist1b) 
    hist2b = hsig2b[index]
    mass =  hist1b.GetName().split('BHBH')[1]
    card = open(os.path.join('datacard_BHBH_M-'+str(mass)+'_shapes_allSysts.txt'), 'w')
    card_content = re.sub('ROOTFILENAME',rootFileName,datacard_template)
    card_content = re.sub('BHBHMBP','BHBH'+str(mass),card_content)
    card_content = re.sub('NOBSCAT1B',str(round(hdata1b.Integral(),3)),card_content)
    card_content = re.sub('NOBSCAT2B',str(round(hdata2b.Integral(),3)),card_content)
    card_content = re.sub('NSIGCAT1B',str(round(hist1b.Integral(),3)),card_content)
    card_content = re.sub('NSIGCAT2B',str(round(hist2b.Integral(),3)),card_content)
    card_content = re.sub('NBKGCAT1B',str(round(hbkg1b.Integral(),3)),card_content)
    card_content = re.sub('NBKGCAT2B',str(round(hbkg2b.Integral(),3)),card_content)
    card.write(card_content)
    card.close() 
    card1b = open(os.path.join('datacard_BHBH_Cat1b_M-'+str(mass)+'_shapes_allSysts.txt'), 'w')
    card1b_content = re.sub('ROOTFILENAME',rootFileName,datacard_template_Cat1b)
    card1b_content = re.sub('BHBHMBP','BHBH'+str(mass),card1b_content)
    card1b_content = re.sub('NOBSCAT1B',str(round(hdata1b.Integral(),3)),card1b_content)
    card1b_content = re.sub('NOBSCAT2B',str(round(hdata2b.Integral(),3)),card1b_content)
    card1b_content = re.sub('NSIGCAT1B',str(round(hist1b.Integral(),3)),card1b_content)
    card1b_content = re.sub('NSIGCAT2B',str(round(hist2b.Integral(),3)),card1b_content)
    card1b_content = re.sub('NBKGCAT1B',str(round(hbkg1b.Integral(),3)),card1b_content)
    card1b_content = re.sub('NBKGCAT2B',str(round(hbkg2b.Integral(),3)),card1b_content)
    card1b.write(card1b_content)
    card1b.close() 
    card2b = open(os.path.join('datacard_BHBH_Cat2b_M-'+str(mass)+'_shapes_allSysts.txt'), 'w')
    card2b_content = re.sub('ROOTFILENAME',rootFileName,datacard_template_Cat2b)
    card2b_content = re.sub('BHBHMBP','BHBH'+str(mass),card2b_content)
    card2b_content = re.sub('NOBSCAT1B',str(round(hdata1b.Integral(),3)),card2b_content)
    card2b_content = re.sub('NOBSCAT2B',str(round(hdata2b.Integral(),3)),card2b_content)
    card2b_content = re.sub('NSIGCAT1B',str(round(hist1b.Integral(),3)),card2b_content)
    card2b_content = re.sub('NSIGCAT2B',str(round(hist2b.Integral(),3)),card2b_content)
    card2b_content = re.sub('NBKGCAT1B',str(round(hbkg1b.Integral(),3)),card2b_content)
    card2b_content = re.sub('NBKGCAT2B',str(round(hbkg2b.Integral(),3)),card2b_content)
    card2b.write(card2b_content)
    card2b.close() 
