#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT

mbps = [500, 600, 700, 800, 1000, 1200, 1500]

tfout = 'ABCDResultTemplate_ABCDana_HT_1.root'

global hist, hist1b, hist2b
tf = ROOT.TFile.Open(tfout, 'update')
tfile = ROOT.TFile.Open('HMassUncert_1sigUp/Final_histograms_BprimebH.root', 'READ')
tfile.cd()
for mass in mbps:
  hist = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_HTAK5')
  hist1b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_1b_HTAK5')
  hist2b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_2b_HTAK5')
  tf.cd()
  hist.Rebin(400)
  hist1b.Rebin(400)
  hist2b.Rebin(400)
  hist.SetName('CatAll_BHBH'+str(mass)+'_MHUp')
  hist1b.SetName('Cat1b_BHBH'+str(mass)+'_MHUp')
  hist2b.SetName('Cat2b_BHBH'+str(mass)+'_MHUp')
  hist.Write()
  hist1b.Write()
  hist2b.Write()
tfile.Close()
tf.Close()

tf = ROOT.TFile.Open(tfout, 'update')
tfile = ROOT.TFile.Open('HMassUncert_1sigDown/Final_histograms_BprimebH.root', 'READ')
tfile.cd()
for mass in mbps:
  hist = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_HTAK5')
  hist1b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_1b_HTAK5')
  hist2b = tfile.Get('BprimeBprimeToBHBHinc_M-'+str(mass)+'__HTSel_2b_HTAK5')
  tf.cd()
  hist.Rebin(400)
  hist1b.Rebin(400)
  hist2b.Rebin(400)
  hist.SetName('CatAll_BHBH'+str(mass)+'_MHDown')
  hist1b.SetName('Cat1b_BHBH'+str(mass)+'_MHDown')
  hist2b.SetName('Cat2b_BHBH'+str(mass)+'_MHDown')
  hist.Write()
  hist1b.Write()
  hist2b.Write()
tfile.Close()
tf.Close()

