#!/bin/tcsh
#if ( $1 == "" ) then
#	echo "Please input the HT cut value"
#	exit
#endif
#source sync.csh
cmsenv
set dataRoots="abcd_JetHT-2012C.root abcd_JetHT-2012D_mis.root abcd_JetHT-2012D.root abcd_JetHT-2012B_mis.root abcd_JetHT-2012B.root abcd_Jet-2012A_mis.root abcd_Jet-2012A.root"	
set sigs="450 500 550 600 650 700 750 800 1000 1200 1500"
#set qcdRoots="abcd_QCD_Pt-300to470.root abcd_QCD_Pt-470to600.root abcd_QCD_Pt-600to800.root abcd_QCD_Pt-800to1000.root abcd_QCD_Pt-1000to1400.root abcd_QCD_Pt-1400to1800.root abcd_QCD_Pt-1800.root"
set qcdRoots="abcd_QCD_Pt-170to300.root abcd_QCD_Pt-300to470.root abcd_QCD_Pt-470to600.root abcd_QCD_Pt-600to800.root abcd_QCD_Pt-800to1000.root abcd_QCD_Pt-1000to1400.root abcd_QCD_Pt-1400to1800.root abcd_QCD_Pt-1800.root"
set ttRoots="abcd_TTJets.root"
#cd /afs/cern.ch/work/j/jtsai/myAna/bpTobH/mywk/CMSSW_5_3_11/src/0501/result/root/abcd/HT$1
#cd /afs/cern.ch/work/j/jtsai/myAna/bpTobH/mywk/CMSSW_5_3_11/src/0501/result/root/abcd
#	foreach sig($sigs)
#		set sigRoot=`echo abcd_BH_$sig.root`
#		rm -f abcd_SumBg_BH_$sig.root 
#		hadd abcd_SumBg_BH_$sig.root $sigRoot $qcdRoots $ttRoots
#	end
	rm -f abcd_OnlySumBg.root
	rm -f abcd_SumQCD.root
	rm -f abcd_AllData.root
	hadd abcd_OnlySumBg.root $qcdRoots $ttRoots
	hadd abcd_SumQCD.root $qcdRoots
	hadd abcd_AllData.root $dataRoots 
#cd -
