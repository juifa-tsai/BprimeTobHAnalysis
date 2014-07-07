#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT

ofname = 'ABCDResultTemplate_ABCDana_HT_02Jul2014.root'
mbps = [500, 600, 700, 800, 900, 1000, 1200, 1500]
global hist, hist1b, hist2b
ifname = sys.argv[1]

dir = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/JobsSig_30June"
fname = "Final_histograms_BprimebH.root"

templ0 = [""]
templUp1 = ["JESUp", "SFbUp", "SFHbUp", "MHUp"]
templUp2 = ["JERUp", "SFlUp", "SFHlUp", "SFCA8Up"]
templUp3 = ["PUUp"]
templDown1 = ["JESDown", "SFbDown", "SFHbDown", "MHDown"]
templDown2 = ["JERDown", "SFlDown", "SFHlDown", "SFCA8Down"]
templDown3 = ["PUDown"]
templ = templ0 + templUp1 + templUp2 + templUp3 + templDown1 + templDown2 + templDown3

fout = ROOT.TFile.Open(ofname, 'RECREATE')

fin = ROOT.TFile.Open(ifname, 'READ')
keys = fin.GetListOfKeys()
key = keys[0]
while key:
  obj = key.ReadObj()
  key = keys.After(key)
  if obj.IsA().InheritsFrom("TH1"):
    if "data_obs" in obj.GetName() or "background" in obj.GetName():
      newobj = obj.Clone()
      fout.cd()
      obj.Write()
fout.Close()

fout = ROOT.TFile.Open(ofname, 'UPDATE')
for itempl in templ: 
  tfile = ROOT.TFile.Open(os.path.join(dir+itempl,fname), 'READ')
  tfile.cd()
  for mass in mbps:
    hist = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_HTAK5')
    hist1b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_1b_HTAK5')
    hist2b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_2b_HTAK5')
    fout.cd()
    histnew = ROOT.TH1D('CatAll_BHBH'+str(mass)+'_'+itempl, '', 9, -250, 2450)
    hist1bnew = ROOT.TH1D('Cat1b_BHBH'+str(mass)+'_'+itempl, '', 9, -250, 2450)
    hist2bnew = ROOT.TH1D('Cat2b_BHBH'+str(mass)+'_'+itempl, '', 9, -250, 2450)
    if itempl == "":
      histnew.SetName(histnew.GetName().strip('_')) 
      hist1bnew.SetName(hist1bnew.GetName().strip('_')) 
      hist2bnew.SetName(hist2bnew.GetName().strip('_')) 
    for bin in range(1, hist.GetNbinsX()+1):
      bincenter = hist.GetBinCenter(bin)
      for newbin in range(1, histnew.GetNbinsX()+1):
        if histnew.GetBinLowEdge(newbin) <= bincenter and bincenter < histnew.GetBinLowEdge(newbin+1):
          hist1bnew.AddBinContent(newbin, hist1b.GetBinContent(bin))
          hist2bnew.AddBinContent(newbin, hist2b.GetBinContent(bin))
    histnew.AddBinContent(8, histnew.GetBinContent(9))
    histnew.AddBinContent(8, histnew.GetBinContent(10))
    histnew.SetBinContent(9, 0)
    histnew.SetBinContent(10, 0)
    hist1bnew.AddBinContent(8, histnew.GetBinContent(9))
    hist1bnew.AddBinContent(8, histnew.GetBinContent(10))
    hist1bnew.SetBinContent(9, 0)
    hist1bnew.SetBinContent(10, 0)
    hist2bnew.AddBinContent(7, histnew.GetBinContent(8))
    hist2bnew.AddBinContent(7, histnew.GetBinContent(9))
    hist2bnew.AddBinContent(7, histnew.GetBinContent(10))
    hist2bnew.SetBinContent(8, 0)
    hist2bnew.SetBinContent(9, 0)
    hist2bnew.SetBinContent(10, 0)
    histnew.Write()
    hist1bnew.Write()
    hist2bnew.Write()
tfile.Close()
fout.Close()
