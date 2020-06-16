#!/usr/bin/env python
import ROOT
import math
import array
import sys
import argparse
import numpy as np
from ROOT import TLorentzVector, AddressOf
from lib2to3 import pytree
import random

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--inputfile', metavar='P', type='string', action='store',
                  default= "",
                  dest='inputFile',
                  help='input file to process')

parser.add_option('--datayear', metavar='P', type='string', action='store',
                  default= "",
                  dest='dataYear',
                  help='input file to process')

parser.add_option('--outputfile', metavar='P', type='string', action='store',
                  default= "",
                  dest='outputFile',
                  help='output file to process')

(options,args) = parser.parse_args()
iFile = options.inputFile
oFile = options.outputFile
year  = str(options.dataYear)
print "year ", year
applyTopPtRew =  False
useDeepCSV = False
#sfFilePath  = "root://cmsxrootd.fnal.gov//store/user/amodak//store/user/amodak/toConndor"
doSys = True
doSkim = False

def Jet_pt(fname, tree, idx, type):
    if ('SingleMuon' in str(fname) or 'SingleElectron' in str(fname) or 'EGamma' in str(fname)):
      return tree.Jet_pt[idx]
    else: 
      if ('nominal' in str(type)):   return tree.Jet_pt_nom[idx]
      elif ('jesUp' in str(type)):   return tree.Jet_pt_jesTotalUp[idx]
      elif ('jesDown' in str(type)): return tree.Jet_pt_jesTotalDown[idx]
      elif ('jerUp' in str(type)):   return tree.Jet_pt_jerUp[idx]
      elif ('jerDown' in str(type)): return tree.Jet_pt_jerDown[idx]

def Jet_mass(fname, tree, idx, type):
    if ('SingleMuon' in str(fname) or 'SingleElectron' in str(fname) or 'EGamma' in str(fname)):
      return tree.Jet_mass[idx]
    else: 
      if ('nominal' in str(type)):   return tree.Jet_mass_nom[idx]
      elif ('jesUp' in str(type)):   return tree.Jet_mass_jesTotalUp[idx]
      elif ('jesDown' in str(type)): return tree.Jet_mass_jesTotalDown[idx]
      elif ('jerUp' in str(type)):   return tree.Jet_mass_jerUp[idx]
      elif ('jerDown' in str(type)): return tree.Jet_mass_jerDown[idx]

def MET_pt(fname, tree, type):
    if ('SingleMuon' in str(fname) or 'SingleElectron' in str(fname) or 'EGamma' in str(fname)):
      return tree.MET_pt
    else: 
      if (year == "2018"): 
        if ('nominal' in str(type)):    return tree.MET_pt_nom
        elif ('jesUp' in str(type)):    return tree.MET_pt_jesTotalUp
        elif ('jesDown' in str(type)):  return tree.MET_pt_jesTotalDown
        elif ('jerUp' in str(type)):    return tree.MET_pt_jerUp
        elif ('jerDown' in str(type)):  return tree.MET_pt_jerDown
      else :
        if ('nominal' in str(type)):    return tree.MET_nom_pt
        elif ('jesUp' in str(type)):    return tree.MET_jesTotalUp_pt
        elif ('jesDown' in str(type)):  return tree.MET_jesTotalDown_pt
        elif ('jerUp' in str(type)):    return tree.MET_jerUp_pt
        elif ('jerDown' in str(type)):  return tree.MET_jerDown_pt

def MET_phi(fname, tree, type):
    if ('SingleMuon' in str(fname) or 'SingleElectron' in str(fname) or 'EGamma' in str(fname)):
      return tree.MET_phi
    else: 
      if (year == "2018"): 
        if ('nominal' in str(type)):    return tree.MET_phi_nom
        elif ('jesUp' in str(type)):    return tree.MET_phi_jesTotalUp
        elif ('jesDown' in str(type)):  return tree.MET_phi_jesTotalDown
        elif ('jerUp' in str(type)):    return tree.MET_phi_jerUp
        elif ('jerDown' in str(type)):  return tree.MET_phi_jerDown
      else :
        if ('nominal' in str(type)):    return tree.MET_nom_phi
        elif ('jesUp' in str(type)):    return tree.MET_jesTotalUp_phi
        elif ('jesDown' in str(type)):  return tree.MET_jesTotalDown_phi
        elif ('jerUp' in str(type)):    return tree.MET_jerUp_phi
        elif ('jerDown' in str(type)):  return tree.MET_jerDown_phi
  
def isMC(fname):
    if ("SingleMuon" in str(fname) or "SingleElectron" in str(fname) or 'EGamma' in str(fname)): 
      return False
    else: 
      return True

def isTTJet(fname):
    if ("TTJets" in str(fname) or "TT_Tune" in str(fname) or "TTTo" in  str(fname)): 
      return True
    else: return False

def btagDeepCSV(btagval, wp):
    flg = 0
    if (year == "2016"): 
        if (wp == "Medium" and btagval > 0.6321): flg = 1
        elif (wp == "Tight" and btagval > 0.8953): flg = 1
    elif (year == "2017"): 
        if (wp == "Medium" and btagval > 0.4941): flg = 1
        elif (wp == "Tight" and btagval > 0.8001): flg = 1
    elif (year == "2018"):
        if (wp == "Medium" and btagval > 0.4184): flg = 1
        elif (wp == "Tight" and btagval > 0.7527): flg = 1
    else: flg = 0
    return flg

def btagDeepFLV(btagval, wp):
    flg = 0
    if (year == "2016"): 
        if (wp == "Medium" and btagval > 0.3093): flg = 1
        elif (wp == "Tight" and btagval > 0.7221): flg = 1
    elif (year == "2017"): 
        if (wp == "Medium" and btagval > 0.3033): flg = 1
        elif (wp == "Tight" and btagval > 0.7489): flg = 1
    elif (year == "2018"):
        if (wp == "Medium" and btagval > 0.2770): flg = 1
        elif (wp == "Tight" and btagval > 0.7264): flg = 1
    else: flg = 0
    return flg

#1: loose (not in 2017), 2:tight, 4: tight lepton veto
def  JetID(idval):
    flag = 0
    if (year == "2016"):
      if (idval >= 2): flag = 1
    elif (year == "2017" or year == "2018"):
      if (idval >= 2): flag = 1
    return flag

def Jet_flavour(entry, idx):
    flavour = -1
    if (abs (entry.Jet_partonFlavour[idx]) ==  5): flavour  = 0
    elif (abs (entry.Jet_partonFlavour[idx]) ==  4): flavour  = 1
    else : flavour  = 2
    return flavour

def Jet_bTagEff(entry, idx):
    eff = 0
    if (abs (entry.Jet_partonFlavour[idx]) ==  5): eff  = 0.80
    elif (abs (entry.Jet_partonFlavour[idx]) ==  4): eff  = 0.20
    else : eff  = 0.01
    return eff

def JetPUID(puid, pt, wp):
    flag = 0
    if (wp == "Loose"): cut = 4
    elif (wp == "Medium"): cut = 6
    elif (wp == "Tight"): cut = 7
    if ((float(puid) >= cut) or (float(pt) > 50)): flag = 1
    return flag

def top_Scaling(ch, nJets):
  #if (year  == "2018" and ch == "Muon"): return (1.0509 - (0.0368*nJets))
  #elif (year  == "2018" and ch == "Electron"): return (1.053 - (0.0421*nJets))
  #else: return 1
  return 1

def st_WJets_Scaling(file, tree):
  weight = 1
  if (isMC(file) and "WJets" in str(file)):
    st = tree.ST_v2
    weight = 1.22 - 0.000219*st
  #return weight
  return 1

def wpt_WJets_Scaling(file, tree):
  weight = 1.0
  if (isMC(file) and "WJets" in str(file)):
    wpt = tree.Lepton_pt + MET_pt(file, tree, "nominal")
    #weight = 1.152 - 0.000321*wpt
    weight = 1.184 - 0.000372*wpt
  #return weight
  return 1

def tm_WJets_Scaling(file, tree, tm):
  weight = 1
  if (isMC(file) and "WJets" in str(file)):
    weight = 1.10 - 0.000197*tm
  #return weight
  return 1

def st_TTJets_Scaling(file, tree):
  weight = 1
  if (isMC(file) and ("TTJets" in str(file) or "TT_Tune" in str(file) or "TTTo" in  str(file))):
    st = tree.ST_v2
    #weight = 1.178 - 0.000298*st
    weight = 1.141 - 0.000231*st
  #return weight
  return 1

def topPt_Reweighting(file, tree, alpha, beta):
    weight = 1
    if (isMC(file)):
      tops = []
      for i in range(0, tree.nGenPart):
        if (abs(tree.GenPart_pdgId[i]) == 6 and abs(tree.GenPart_status[i]) == 22):
          top_4vec = ROOT.TLorentzVector()
          top_4vec.SetPtEtaPhiM(tree.GenPart_pt[i], tree.GenPart_eta[i], tree.GenPart_phi[i], tree.GenPart_mass[i])
          tops.append(top_4vec)
      if (len(tops) == 2):
         pt1 = tops[0].Pt()
         pt2 = tops[1].Pt()
         weight = 0.89*1.34*math.sqrt(math.exp(alpha - beta * pt1) * math.exp(alpha - beta * pt2))
    return weight

def W_GeneratorP4_weight(file, tree, alpha, beta):
    weight = 1
    if (isMC(file)):
      for i in range(0, tree.nGenPart):
        if (abs(tree.GenPart_pdgId[i]) == 24 and tree.GenPart_status[i] == 22):
          w_4vec = ROOT.TLorentzVector()
          w_4vec.SetPtEtaPhiM(tree.GenPart_pt[i], tree.GenPart_eta[i], tree.GenPart_phi[i], tree.GenPart_mass[i])
          #new_momentum =  beta*w_4vec.P()
          #new_pt = new_momentum/math.cosh(tree.GenPart_eta[i])
          weight = alpha + beta*w_4vec.P()
    return weight

def applyHEM(file, entry, j):
  flag = 1
  if (year == "2018"):
    if ((entry.Jet_eta[j] < -1.3 and entry.Jet_eta[j] > -3.0 and entry.Jet_phi[j] > -1.57 and entry.Jet_phi[j]  < -0.87) or (entry.Lepton_eta < -1.3 and entry.Lepton_eta > -3.0 and entry.Lepton_phi > -1.57 and entry.Lepton_phi  < -0.87)): foundLJ = 1
    else: foundLJ = 0
    #Found Lepton or Jet in HEM 15/16
    if (foundLJ):
      #print "Jet phi ", entry.Jet_phi[j], " Jet eta ", entry.Jet_eta[j], " Lep phi ", entry.Lepton_phi, " Lep eta", entry.Lepton_eta
      if (isMC(file)):
        #Affected lumi 39/60 = 0.65
        #ran_proportion = random.uniform(0, 1)
        #print "ran_proportion", ran_proportion
        #if (ran_proportion < 0.65):
        flag = 0
      else:
        if ("SingleMuon_2018C" in str(file) or "SingleMuon_2018D" in str(file) or "EGamma_2018C" in str(file) or "EGamma_2018D" in str(file)):
          flag = 0
  return flag

def trig_SF(file, tree, pt, eta, channel):
  sf = 1
  if (year == "2018" and channel  == "Muon" and (isMC(file))):
    sf = tree.GetBinContent(tree.GetXaxis().FindBin(pt), tree.GetYaxis().FindBin(abs(eta)))
  return sf

def MET_Filters(file, tree):
    if (year == "2016"):
      if (isMC(file)):
        if (tree.Flag_goodVertices > 0 and tree.Flag_globalSuperTightHalo2016Filter > 0 and tree.Flag_HBHENoiseFilter > 0 and  tree.Flag_HBHENoiseIsoFilter > 0 and tree.Flag_EcalDeadCellTriggerPrimitiveFilter > 0 and tree.Flag_BadPFMuonFilter > 0): 
          return 1
        else:
          return 0
      else:
        if (tree.Flag_goodVertices > 0 and tree.Flag_globalSuperTightHalo2016Filter > 0 and tree.Flag_HBHENoiseFilter > 0 and  tree.Flag_HBHENoiseIsoFilter > 0 and tree.Flag_EcalDeadCellTriggerPrimitiveFilter > 0 and tree.Flag_BadPFMuonFilter > 0 and tree.Flag_eeBadScFilter > 0): 
          return 1
        else:
          return 0
    elif (year == "2017" or year == "2018"):
      if (isMC(file)):
        if (tree.Flag_goodVertices > 0 and tree.Flag_globalSuperTightHalo2016Filter > 0 and tree.Flag_HBHENoiseFilter > 0 and  tree.Flag_HBHENoiseIsoFilter > 0 and tree.Flag_EcalDeadCellTriggerPrimitiveFilter > 0 and tree.Flag_BadPFMuonFilter > 0): 
          return 1
        else:
          return 0
      else:
        if (tree.Flag_goodVertices > 0 and tree.Flag_globalSuperTightHalo2016Filter > 0 and tree.Flag_HBHENoiseFilter > 0 and  tree.Flag_HBHENoiseIsoFilter > 0 and tree.Flag_EcalDeadCellTriggerPrimitiveFilter > 0 and tree.Flag_BadPFMuonFilter > 0 and tree.Flag_eeBadScFilter > 0): 
          return 1
        else:
          return 0
    else: return 0

def applyBTag (isBTagged, sf, eff, ran):
    btagFlag = isBTagged
    if (sf == 1): return btagFlag
    if (sf > 1):
      if not (isBTagged == 1):
        #fraction of jets that need to be upgraded
        mistagPercent = (1.0 - sf)/(1.0 - (1.0/eff))
        #upgrade to tagged
        if (ran < mistagPercent): btagFlag = 1 
    else:
      #downgrade tagged to untagged
      if (isBTagged == 1 and ran > sf): btagFlag = 0
    return btagFlag

def applyHighPtMuonSF (file, entry, pt, eta):
    sf = 1
    data_eff = 1
    mc_eff  = 1
    if (abs (eta) < 1.6 and pt < 100):
      mc_eff = 0.989 - 2.399*(10.0**-6)*pt 
    if (abs (eta) < 1.6 and pt > 100):
      data_eff = 0.9828 - 1.947*(10.0**-5)*pt 
    if (abs (eta) > 1.6 and abs (eta) < 2.4 and pt < 275):
      mc_eff = 0.9974 - 1.721*(10.0**-5)*pt 
    if (abs (eta) > 1.6 and abs (eta) < 2.4 and pt > 275):
      data_eff = 0.9893 - 3.666*(10.0**-5)*pt 
    sf = float(data_eff/mc_eff)
    return sf

def cutFlow(file, entry, crType, leadjetBTag, HT_Central, nBTag, MassT, MET_pt, LeadJet_pt):
  if (crType == "WJets"):
    if not (entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1): return 1
    if not (MET_pt > 60): return 2
    if not (MET_Filters(file, entry) > 0): return 3
    if not (LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4): return 4
    if not (leadjetBTag < 1): return 5
    if not (entry.DR_LepleadJet > 0.5): return 6
    if not (entry.DR_LepClosestJet > -1): return 7
    if not (MassT > 40): return 8 ##Changed 40->50
    if not (nBTag == 0): return 9
    if not (entry.ST_v2 > 500): return 10
    if not (HT_Central > 500): return 11
    if not (abs(entry.DPHI_LepMet) < 0.5): return 12
    if not (abs(entry.DPHI_LepleadJet) > 2.0): return 13
    if not (entry.bVsW_ratio < 999): return 14
    if not ((entry.nVetoMuons + entry.nVetoElectrons) == 1): return 15
    if not (entry.FwdJetPt > -999 and abs(entry.FwdJetEta) > -999): return 16  ##changed
    return 17
  if (crType == "Signal"):
    if not (entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1): return 1
    if not (MET_pt > 60): return 2
    if not (MET_Filters(file, entry) > 0): return 3
    if not (LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4): return 4
    if not (leadjetBTag >= 1): return 5
    if not (entry.DR_LepleadJet > 0.5): return 6
    if not (entry.DR_LepClosestJet > 1.5): return 7
    if not (MassT > 40): return 8
    if not (nBTag >= 1): return 9
    if not (entry.ST_v2 > 700): return 10
    if not (HT_Central > 500): return 11
    if not (abs(entry.DPHI_LepMet) < 0.5): return 12
    if not (abs(entry.DPHI_LepleadJet) > 2.0): return 13
    if not (entry.bVsW_ratio < 1.4): return 14
    if not ((entry.nVetoMuons + entry.nVetoElectrons) == 1): return 15
    if not (entry.FwdJetPt > 30 and abs(entry.FwdJetEta) > 2.4): return 16 ##Changed
    return 17
  if (crType == "TTJets"):
    if not (entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1): return 1
    if not (MET_pt > 60): return 2
    if not (MET_Filters(file, entry) > 0): return 3
    if not (LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4): return 4
    if not (leadjetBTag >= 1): return 5
    if not (entry.DR_LepleadJet > -1): return 6
    if not (entry.DR_LepClosestJet > 0.5 and entry.DR_LepClosestJet < 1.5): return 7
    if not (MassT > 0): return 8
    if not (nBTag >= 2): return 9
    if not (entry.ST_v2 > 300): return 10
    if not (HT_Central > 500): return 11
    if not (abs(entry.DPHI_LepMet) < 999): return 12
    if not (abs(entry.DPHI_LepleadJet) > -999): return 13
    if not (entry.bVsW_ratio < 999): return 14
    if not ((entry.nVetoMuons + entry.nVetoElectrons) == 1): return 15
    if not (entry.FwdJetPt > -99 and abs(entry.FwdJetEta) > -99): return 16  ##Changed
    return 17
  if (crType == "MultiJets"):
    if not (entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1): return 1
    if not (MET_pt > 60): return 2
    if not (MET_Filters(file, entry) > 0): return 3
    if not (LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4): return 4
    if not (leadjetBTag >= 1): return 5
    if not (entry.DR_LepleadJet > -1): return 6
    if not (entry.DR_LepClosestJet > -1.0): return 7
    if not (MassT < 9999): return 8
    if not (nBTag >= 1): return 9
    if not (entry.ST_v2 > 300): return 10
    if not (HT_Central > 200): return 11
    if not (abs(entry.DPHI_LepMet) < 999): return 12
    if not (abs(entry.DPHI_LepleadJet) > -999): return 13
    if not (entry.bVsW_ratio < 999): return 14
    if not ((entry.nVetoMuons + entry.nVetoElectrons) == 1): return 15
    #if not (entry.FwdJetPt > -99 and abs(entry.FwdJetEta) > -99): return 16
    if not (entry.nCentralJets >= 4): return 16
    return 17
  if (crType == "PreSig"):
    if not (entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1): return 1
    if not (MET_pt > 60): return 2
    if not (MET_Filters(file, entry) > 0): return 3
    if not (LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4): return 4
    if not (leadjetBTag >= 1): return 5
    if not (entry.DR_LepleadJet > 0.5): return 6
    if not (entry.DR_LepClosestJet > -1): return 7
    if not (MassT > 40): return 8
    if not (nBTag >= 1): return 9
    if not (entry.ST_v2 > 500): return 10
    if not (HT_Central > 500): return 11
    if not (abs(entry.DPHI_LepMet) < 999): return 12
    if not (abs(entry.DPHI_LepleadJet) > -999): return 13
    if not (entry.bVsW_ratio < 999): return 14
    if not ((entry.nVetoMuons + entry.nVetoElectrons) == 1): return 15
    if not (entry.FwdJetPt > -99 and abs(entry.FwdJetEta) > -99): return 16
    return 17

def fillHisto(cntR, channeL, file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt, MET_phi, LeadJet_pt, sys, evtwt):
    iso03  = 999
    iso04  = 999
    if (channeL == "Mu"):
      iso03  = entry.Muon_pfRelIso03_all[entry.Lepton_idx]
      iso04  = entry.Muon_pfRelIso04_all[entry.Lepton_idx]
    elif (channeL == "Ele"):
      iso03  = entry.Electron_pfRelIso03_all[entry.Lepton_idx]
      
    hmap[str(channeL)+"_"+str(cntR)+"_LepPt_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_pt, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_LepPhi_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_phi, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_LepEta_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_eta, evtwt)
    if (HEMwt == 0):
      hmap[str(channeL)+"_"+str(cntR)+"_LepPhi_HEMf_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_phi, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_LepEta_HEMf_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_eta, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_JetPhi_HEMf_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_phi, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_JetEta_HEMf_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_eta, evtwt)
    elif (HEMwt == 1):
      hmap[str(channeL)+"_"+str(cntR)+"_LepPhi_HEMp_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_phi, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_LepEta_HEMp_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_eta, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_JetPhi_HEMp_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_phi, evtwt)
      hmap[str(channeL)+"_"+str(cntR)+"_JetEta_HEMp_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_eta, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_MET_"+str(sys)+"_"+str(itr)].Fill(MET_pt, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_METphi_"+str(sys)+"_"+str(itr)].Fill(MET_phi, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_ST_v2_"+str(sys)+"_"+str(itr)].Fill(entry.ST_v2, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_TranMom_"+str(sys)+"_"+str(itr)].Fill((entry.ST_v2+HT_Central-LeadJet_pt), evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_ST_"+str(sys)+"_"+str(itr)].Fill(entry.ST, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_WPt_"+str(sys)+"_"+str(itr)].Fill(entry.Lepton_pt + MET_pt, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_MT_"+str(sys)+"_"+str(itr)].Fill(MassT, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_HT_"+str(sys)+"_"+str(itr)].Fill(HT_Central, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_DPHI_"+str(sys)+"_"+str(itr)].Fill(entry.DPHI_LepleadJet, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_DPHILepMet_"+str(sys)+"_"+str(itr)].Fill(entry.DPHI_LepMet, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_DR_"+str(sys)+"_"+str(itr)].Fill(entry.DR_LepClosestJet, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_DR_LepleadJet_"+str(sys)+"_"+str(itr)].Fill(entry.DR_LepleadJet, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_JetPt_"+str(sys)+"_"+str(itr)].Fill(LeadJet_pt, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_JetPhi_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_phi, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_JetEta_"+str(sys)+"_"+str(itr)].Fill(entry.LeadJet_eta, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_NBTags_DeepCSV_"+str(sys)+"_"+str(itr)].Fill(entry.nBTagMed_DeepCSV, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_NBTags_DeepFLV_"+str(sys)+"_"+str(itr)].Fill(nBTag, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_NJets_"+str(sys)+"_"+str(itr)].Fill(entry.nJets, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_NCentralJets_"+str(sys)+"_"+str(itr)].Fill(nCentralJets_v2, evtwt)
    #hmap[str(channeL)+"_"+str(cntR)+"_NCentralJets_"+str(sys)+"_"+str(itr)].Fill(entry.nCentralJets, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_NFwdJets_"+str(sys)+"_"+str(itr)].Fill(entry.nFwdJets, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_DPHIMetJet_"+str(sys)+"_"+str(itr)].Fill(entry.DPHI_JetMet, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_FwdJetEta_"+str(sys)+"_"+str(itr)].Fill(entry.FwdJetEta, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_FwdJetPt_"+str(sys)+"_"+str(itr)].Fill(entry.FwdJetPt, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_FwdJetPt_L_"+str(sys)+"_"+str(itr)].Fill(FwdJetPt_L, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_FwdJetPt_M_"+str(sys)+"_"+str(itr)].Fill(FwdJetPt_M, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_FwdJetPt_T_"+str(sys)+"_"+str(itr)].Fill(FwdJetPt_T, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_Mass_"+str(sys)+"_"+str(itr)].Fill(entry.Mass, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_Mass_v2_"+str(sys)+"_"+str(itr)].Fill(entry.Mass_v2, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_RelIso03_"+str(sys)+"_"+str(itr)].Fill(iso03, evtwt)
    hmap[str(channeL)+"_"+str(cntR)+"_RelIso04_"+str(sys)+"_"+str(itr)].Fill(iso04, evtwt)
    if (nBTag == 1):
       hmap[str(channeL)+"_"+str(cntR)+"_Mass1b_"+str(sys)+"_"+str(itr)].Fill(entry.Mass, evtwt)
       hmap[str(channeL)+"_"+str(cntR)+"_Mass1b_v2_"+str(sys)+"_"+str(itr)].Fill(entry.Mass_v2, evtwt)
    elif (nBTag > 1):
       hmap[str(channeL)+"_"+str(cntR)+"_Mass2b_"+str(sys)+"_"+str(itr)].Fill(entry.Mass, evtwt)
       hmap[str(channeL)+"_"+str(cntR)+"_Mass2b_v2_"+str(sys)+"_"+str(itr)].Fill(entry.Mass_v2, evtwt)


for i in  range(1):
  file = ROOT.TFile.Open(iFile)
  print "Opened EOS file to process ..."
  tree = file.Get("Events")
  #canv = ROOT.TCanvas()
  outfile = oFile
  outputfile = ROOT.TFile(outfile, 'RECREATE')
  #Skim Tree
  if (doSkim):
    sktree = ROOT.TTree( 'Skim', 'Skim Tree' )
    Event_flag = array.array('i', [0])
    Lepton_pt = array.array('f', [0])
    Lepton_eta = array.array('f', [0])
    Lepton_phi = array.array('f', [0])
    Lepton_mass = array.array('f', [0]) 
    LeadJet_pt = array.array('f', [0]) 
    LeadJet_eta = array.array('f', [0])
    LeadJet_phi = array.array('f', [0])
    LeadJet_mass = array.array('f', [0])
    LeadJet_btagDeepCSV = array.array('f', [0])  
    LeadJet_btagDeepFLV = array.array('f', [0])  
    LeadJet_btag = array.array('f', [0])  
    MET = array.array('f', [0])
    MT = array.array('f', [0])
    HT = array.array('f', [0])
    ST = array.array('f', [0])
    ST_v2 = array.array('f', [0])
    nVetoElectrons = array.array('i', [0])
    nVetoMuons = array.array('i', [0])
    Mass = array.array('f', [0])
    Mass_v2 = array.array('f', [0])
    DR_LepClosestJet = array.array('f', [0])
    DR_LepleadJet = array.array('f', [0])
    DPHI_LepleadJet = array.array('f', [0])
    DPHI_LepMet = array.array('f', [0])
    DPHI_JetMet = array.array('f', [0])
    Angle_MuZ_Jet = array.array('f', [0])
    Angle_MuJet_Met = array.array('f', [0])
    bVsW_ratio = array.array('f', [0])
    nJets = array.array('i', [0])
    nCentralJets = array.array('i', [0])
    nFwdJets = array.array('i', [0])
    nBTagMed_DeepCSV = array.array('i', [0])
    nBTagTig_DeepCSV = array.array('i', [0])
    nBTagMed_DeepFLV = array.array('i', [0])
    nBTagTig_DeepFLV = array.array('i', [0])
    FwdJetPt = array.array('f', [0])
    FwdJetEta = array.array('f', [0])
    FwdJetPhi = array.array('f', [0])
    STScale_TT = array.array('f', [0])
    WPTScale_WJ = array.array('f', [0])
    EvtWt = array.array('f', [0])

    sktree.Branch("Event_flag", Event_flag,'Event_flag/I')
    sktree.Branch("Lepton_pt", Lepton_pt, 'Lepton_pt/F')
    sktree.Branch("Lepton_eta", Lepton_eta, 'Lepton_eta/F')
    sktree.Branch("Lepton_phi", Lepton_phi, 'Lepton_phi/F')
    sktree.Branch("Lepton_mass", Lepton_mass, 'Lepton_mass/F')
    sktree.Branch("LeadJet_pt", LeadJet_pt, 'LeadJet_pt/F')
    sktree.Branch("LeadJet_eta", LeadJet_eta, 'LeadJet_eta/F')
    sktree.Branch("LeadJet_phi", LeadJet_phi, 'LeadJet_phi/F')
    sktree.Branch("LeadJet_mass", LeadJet_mass, 'LeadJet_mass/F')
    sktree.Branch("LeadJet_btagDeepCSV", LeadJet_btagDeepCSV, 'LeadJet_btagDeepCSV/F')
    sktree.Branch("LeadJet_btagDeepFLV", LeadJet_btagDeepFLV, 'LeadJet_btagDeepFLV/F')
    sktree.Branch("LeadJet_btag", LeadJet_btag, 'LeadJet_btag/F')
    sktree.Branch("MET", MET, 'MET/F')
    sktree.Branch("MT", MT, 'MT/F')
    sktree.Branch("HT", HT, 'HT/F')
    sktree.Branch("ST", ST, 'ST/F')
    sktree.Branch("ST_v2", ST_v2, 'ST_v2/F')
    sktree.Branch("nVetoElectrons", nVetoElectrons, 'nVetoElectrons/I')
    sktree.Branch("nVetoMuons", nVetoMuons, 'nVetoMuons/I')
    sktree.Branch("Mass", Mass, 'Mass/F')
    sktree.Branch("Mass_v2", Mass_v2, 'Mass_v2/F')
    sktree.Branch("DR_LepClosestJet", DR_LepClosestJet, 'DR_LepClosestJet/F')
    sktree.Branch("DR_LepleadJet", DR_LepleadJet, 'DR_LepleadJet/F')
    sktree.Branch("DPHI_LepleadJet", DPHI_LepleadJet, 'DPHI_LepleadJet/F')
    sktree.Branch("DPHI_LepMet", DPHI_LepMet, 'DPHI_LepMet/F')
    sktree.Branch("DPHI_JetMet", DPHI_JetMet, 'DPHI_JetMet/F')
    sktree.Branch("Angle_MuZ_Jet", Angle_MuZ_Jet, 'Angle_MuZ_Jet/F')
    sktree.Branch("Angle_MuJet_Met", Angle_MuJet_Met, 'Angle_MuJet_Met/F')
    sktree.Branch("bVsW_ratio", bVsW_ratio, 'bVsW_ratio/F')
    sktree.Branch("nJets", nJets, 'nJets/I')
    sktree.Branch("nCentralJets", nCentralJets, 'nCentralJets/I')
    sktree.Branch("nFwdJets", nFwdJets, 'nFwdJets/I')
    sktree.Branch("nBTagMed_DeepCSV", nBTagMed_DeepCSV, 'nBTagMed_DeepCSV/I')
    sktree.Branch("nBTagTig_DeepCSV", nBTagTig_DeepCSV, 'nBTagTig_DeepCSV/I')
    sktree.Branch("nBTagMed_DeepFLV", nBTagMed_DeepFLV, 'nBTagMed_DeepFLV/I')
    sktree.Branch("nBTagTig_DeepFLV", nBTagTig_DeepFLV, 'nBTagTig_DeepFLV/I')
    sktree.Branch("FwdJetPt", FwdJetPt, 'FwdJetPt/F')
    sktree.Branch("FwdJetEta", FwdJetEta, 'FwdJetEta/F')
    sktree.Branch("FwdJetPhi", FwdJetPhi, 'FwdJetPhi/F')
    sktree.Branch("STScale_TT", STScale_TT, 'EvtWt/F')
    sktree.Branch("WPTScale_WJ", WPTScale_WJ, 'EvtWt/F')
    sktree.Branch("EvtWt", EvtWt, 'EvtWt/F')


  #Histograms
  Entry           = ROOT.TH1D("Entry", "", 1, 0.5, 1.5)
  Entry_nowt      = ROOT.TH1D("Entry_nowt", "", 1, 0.5, 1.5)
  Entry_genwt      = ROOT.TH1D("Entry_genwt", "", 1, 0.5, 1.5)
  effmap_pt_eta_b_all = ROOT.TH2D("effmap_pt_eta_b_all", "", 100, 30, 1030, 10, -3, 3)
  effmap_pt_eta_c_all = ROOT.TH2D("effmap_pt_eta_c_all", "", 100, 30, 1030, 10, -3, 3)
  effmap_pt_eta_l_all = ROOT.TH2D("effmap_pt_eta_l_all", "", 100, 30, 1030, 10, -3, 3)
  effmap_pt_eta_b_btagged = ROOT.TH2D("effmap_pt_eta_b_btagged", "", 100, 30, 1030, 10, -3, 3)
  effmap_pt_eta_c_btagged = ROOT.TH2D("effmap_pt_eta_c_btagged", "", 100, 30, 1030, 10, -3, 3)
  effmap_pt_eta_l_btagged = ROOT.TH2D("effmap_pt_eta_l_btagged", "", 100, 30, 1030, 10, -3, 3)
  Mu_TTJets_Mass_lead_sublead_bin0 = ROOT.TH2D("Mu_TTJets_Mass_lead_sublead_bin0", "", 1000, 0, 1000, 1000, 0, 1000)
  Mu_TTJets_Mass_lead_sublead_bin1 = ROOT.TH2D("Mu_TTJets_Mass_lead_sublead_bin1", "", 1000, 0, 1000, 1000, 0, 1000)
  Mu_TTJets_Mass_lead_sublead_bin2 = ROOT.TH2D("Mu_TTJets_Mass_lead_sublead_bin2", "", 1000, 0, 1000, 1000, 0, 1000)
  Mu_TTJets_Mass_lead_sublead_bin3 = ROOT.TH2D("Mu_TTJets_Mass_lead_sublead_bin3", "", 1000, 0, 1000, 1000, 0, 1000)
  Ele_TTJets_Mass_lead_sublead_bin0 = ROOT.TH2D("Ele_TTJets_Mass_lead_sublead_bin0", "", 1000, 0, 1000, 1000, 0, 1000)
  Ele_TTJets_Mass_lead_sublead_bin1 = ROOT.TH2D("Ele_TTJets_Mass_lead_sublead_bin1", "", 1000, 0, 1000, 1000, 0, 1000)
  Ele_TTJets_Mass_lead_sublead_bin2 = ROOT.TH2D("Ele_TTJets_Mass_lead_sublead_bin2", "", 1000, 0, 1000, 1000, 0, 1000)
  Ele_TTJets_Mass_lead_sublead_bin3 = ROOT.TH2D("Ele_TTJets_Mass_lead_sublead_bin3", "", 1000, 0, 1000, 1000, 0, 1000)
  Mu_WJets_HT_vs_ST = ROOT.TH2D("Mu_WJets_HT_vs_ST", "", 300, 0, 3000, 300, 0, 3000)
  Ele_WJets_HT_vs_ST = ROOT.TH2D("Ele_WJets_HT_vs_ST", "", 300, 0, 3000, 300, 0, 3000)
  Jet_eta_phi_2D = ROOT.TH2D("Jet_eta_phi_2D", "", 50, -2.5, 2.5, 70, -3.5, 3.5)
  Lep_eta_phi_2D = ROOT.TH2D("Lep_eta_phi_2D", "", 50, -2.5, 2.5, 70, -3.5, 3.5)

  hmap = {}
  channels = ["Mu", "Ele"]
  controlR = ["WJets", "TTJets", "Signal", "PreSig"]
  #controlR = ["wjet", "ttjet", "signal", "presig", "multijet"]
  variables = {"LepPt":  [100, 0, 1000], "LepPhi":  [100, -5, 5], "LepPhi_HEMf":  [100, -5, 5], "LepPhi_HEMp":  [100, -5, 5], "MET": [100, 0, 1000], "METphi": [100, -5, 5], "LepEta": [120, -3, 3], "LepEta_HEMf": [120, -3, 3], "LepEta_HEMp": [120, -3, 3], "ST": [300, 0, 3000], "ST_v2": [300, 0, 3000], "MT": [200, 0, 200], "HT": [300, 0, 3000], "DPHI": [100, -5, 5], "DPHILepMet": [100, -5, 5], "DR": [100, 0, 5], "DR_LepleadJet": [100, 0, 5], "JetPt": [200, 0, 2000], "JetPhi": [100, -5, 5], "JetEta": [100, -5, 5], "JetPhi_HEMp": [100, -5, 5], "JetEta_HEMp": [100, -5, 5], "JetPhi_HEMf": [100, -5, 5], "JetEta_HEMf": [100, -5, 5], "NBTags_DeepCSV": [5, -0.5, 4.5], "NBTags_DeepFLV": [5, -0.5, 4.5], "NJets": [20, 0, 20], "NCentralJets": [20, 0, 20], "NFwdJets": [20, 0, 20], "DPHIMetJet": [100, -5, 5], "FwdJetEta": [100, -5, 5], "FwdJetPt": [500, 0, 500], "Mass": [300, 0, 3000], "Mass_v2": [300, 0, 3000], "Mass1b": [300, 0, 3000], "Mass1b_v2": [300, 0, 3000], "Mass2b": [300, 0, 3000], "Mass2b_v2": [300, 0, 3000], "RelIso03": [100, 0, 0.5], "RelIso04": [100, 0, 0.5], "WPt": [200, 0, 2000], "TranMom": [600, 0, 6000], "FwdJetPt_L": [500, 0, 500], "FwdJetPt_M": [500, 0, 500], "FwdJetPt_T": [500, 0, 500]}

  variables_top = {"Pt_lmj_select": [100, 0, 1000], "Pt_lmj_select2": [100, 0, 1000], "Pt_lmj_select_alpha_up": [100, 0, 1000], "Pt_lmj_select_alpha_down": [100, 0, 1000], "Pt_lmj_select_beta_up": [100, 0, 1000], "Pt_lmj_select_beta_down": [100, 0, 1000], "Top_Score": [4, -0.5, 3.5], "Pt_lmj_select_beta_up2": [100, 0, 1000], "Pt_lmj_select_beta_down2": [100, 0, 1000], "ST_v2_def2": [300, 0, 3000], "ST_v2_beta_up2": [300, 0, 3000], "ST_v2_beta_down2": [300, 0, 3000]}

  #levels = [4, 11, 16]
  levels = [16]
  systematicS = ["nominal"]
  if (doSys):
    systematicS = ["nominal", "PileupUp", "PileupDown", "BTagSFUp", "BTagSFDown", "topptweightUp", "topptweightDown", "jerUp", "jerDown", "jesUp", "jesDown", "LHEScale1", "LHEScale2", "LHEScale3", "LHEScale4", "LHEScale6", "LHEScale8", "LHEScaleUpWeight", "LHEScaleDownWeight", "TTJets_ST_Scaling"]
    
  for ch in range(len(channels)):
    for cr in range(len(controlR)):
      hname = channels[ch]+"_"+controlR[cr]+"_Counter"
      hnamestr = hname
      hname = ROOT.TH1D(str(hnamestr), "", 21, -0.5, 20.5)
      hmap.update({str(hnamestr) : hname})
      for key in variables_top:
        hname = channels[ch]+"_"+controlR[cr]+"_"+key
        hnamestr = hname
        hname = ROOT.TH1D(str(hnamestr), "", variables_top[key][0], variables_top[key][1], variables_top[key][2])
        hmap.update({str(hnamestr) : hname})
      for hist in levels:
        for key in variables:
          for item in systematicS:
            hname = channels[ch]+"_"+controlR[cr]+"_"+key+"_"+str(item)+"_"+str(hist)
            hnamestr = hname
            hname = ROOT.TH1D(str(hnamestr), "", variables[key][0], variables[key][1], variables[key][2])
            hmap.update({str(hnamestr) : hname})

  if (isMC(file)): 
    histo = file.Get("autoPU")
    Entry_nowt.Fill(1, float(histo.GetEntries()))
    runs = file.Get("Runs")
    for itr in runs:
      if (year == "2018"):
        Entry_genwt.Fill(1, itr.genEventSumw_)
      else :
        Entry_genwt.Fill(1, itr.genEventSumw)
    #Open SF root files here
    if (year == "2018"): 
      file_muTrigBefore = ROOT.TFile.Open("root://cmsxrootd.fnal.gov//store/user/amodak/toConndor/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root")
      tree_muTrigBefore = file_muTrigBefore.Get("IsoMu24_PtEtaBins/pt_abseta_ratio")
      file_muTrigAfter = ROOT.TFile.Open("root://cmsxrootd.fnal.gov//store/user/amodak/toConndor/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root")
      tree_muTrigAfter = file_muTrigAfter.Get("IsoMu24_PtEtaBins/pt_abseta_ratio")

  #Event loop starts
  for entry in tree:
    if (isMC(file)): 
      Entry.Fill(1, float(entry.genWeight)*float(entry.puWeight))
      #runs = file.Get("Runs")
      #for itr in runs:
      #  if (year == "2018"):
      #    for nScale in range(itr.nLHEScaleSumw_):
      #      print "Index: ", nScale, ", nLHEScale: ", itr.LHEScaleSumw_[nScale]
    else: 
      Entry.Fill(1, 1)
    MassT = -1
    Lep4vec, Met4vec = TLorentzVector(), TLorentzVector()
    Lep4vec.SetPtEtaPhiM(entry.Lepton_pt, entry.Lepton_eta, entry.Lepton_phi, entry.Lepton_mass)
    Met4vec.SetPtEtaPhiM(MET_pt(file, entry, "nominal"), 0, MET_phi(file, entry, "nominal"), 0)
    MassT = math.sqrt(2*entry.Lepton_pt*MET_pt(file, entry, "nominal")*(1 - math.cos(Lep4vec.DeltaPhi(Met4vec))))
 
    leadjetBTag = -999
    nBTag  = 0
    #Tuned to "Medium" WP, change here if needed
    if (useDeepCSV): 
      nBTag = entry.nBTagMed_DeepCSV
      leadjetBTag = btagDeepCSV(entry.LeadJet_btag_DeepCSV, "Medium")
    else:
      nBTag = entry.nBTagMed_DeepFLV
      leadjetBTag = btagDeepFLV(entry.LeadJet_btag_DeepFLV, "Medium")

    HT_Central =  0
    nCentralJets_v2 = 0
    dR_lepClosestJet_v2 = 999
    lepton4vec = TLorentzVector()
    lepton4vec.SetPtEtaPhiM(entry.Lepton_pt, entry.Lepton_eta, entry.Lepton_phi, entry.Lepton_mass)
    HEMwt  = 1
    FwdJetPt_L = -999
    FwdJetPt_M = -999
    FwdJetPt_T = -999
    fwdJeteta_L =  0
    fwdJeteta_M =  0
    fwdJeteta_T =  0

    for j in range(0, entry.nJet):
       if (Jet_pt(file, entry, j, "nominal") > 30 and JetID(entry.Jet_jetId[j]) >= 1): ##Changed
         #Check HEM for all central jets
         if (HEMwt == 1 and applyHEM(file, entry, j) == 0): HEMwt = 0
         jet4vec = TLorentzVector()
         jet4vec.SetPtEtaPhiM(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j], entry.Jet_phi[j], Jet_mass(file, entry, j, "nominal"))
         if (jet4vec.DeltaR(Lep4vec) > 0.5 and abs(entry.Jet_eta[j]) < 2.4 and JetPUID(entry.Jet_puId[j], Jet_pt(file, entry, j, "nominal"), "Tight") >= 1):
           HT_Central += jet4vec.Pt()
           nCentralJets_v2 += 1

         if (jet4vec.DeltaR(lepton4vec) < dR_lepClosestJet_v2 and Jet_pt(file, entry, j, "nominal") > 40 and JetPUID(entry.Jet_puId[j], Jet_pt(file, entry, j, "nominal"), "Tight") >= 1 and abs(entry.Jet_eta[j]) < 2.4): 
           dR_lepClosestJet_v2 = jet4vec.DeltaR(lepton4vec)
         
         if (abs(entry.Jet_eta[j]) > abs(fwdJeteta_L)  and  abs(entry.Jet_eta[j]) > 2.4 and JetPUID(entry.Jet_puId[j], Jet_pt(file, entry, j, "nominal"), "Loose") >= 1):
           fwdJeteta_L = entry.Jet_eta[j]
           FwdJetPt_L = Jet_pt(file, entry, j, "nominal")
         if (abs(entry.Jet_eta[j]) > abs(fwdJeteta_L)  and  abs(entry.Jet_eta[j]) > 2.4 and JetPUID(entry.Jet_puId[j], Jet_pt(file, entry, j, "nominal"), "Medium") >= 1):
           fwdJeteta_M = entry.Jet_eta[j]
           FwdJetPt_M = Jet_pt(file, entry, j, "nominal")
         if (abs(entry.Jet_eta[j]) > abs(fwdJeteta_L)  and  abs(entry.Jet_eta[j]) > 2.4 and JetPUID(entry.Jet_puId[j], Jet_pt(file, entry, j, "nominal"), "Tight") >= 1):
           fwdJeteta_T = entry.Jet_eta[j]
           FwdJetPt_T = Jet_pt(file, entry, j, "nominal")

    btagSF_0tag = 1
    btagSF_1tag = 1
    btagSF_2tag = 1
    btagSF = []
    btagSF_evt = 1
    btagSF_evt_tree = 1
    btagSF_evt_treeUp = 1
    btagSF_evt_treeDown = 1

    #Event weights
    evtwt = 1
    evtwt_nominal = 1
    evtwt_PileupUp = 1
    evtwt_PileupDown = 1
    evtwt_BTagSFUp = 1
    evtwt_BTagSFDown = 1
    evtwt_topptweightUp = 1
    evtwt_topptweightDown = 1
    evtwt_lhescale1 = 1   
    evtwt_lhescale2 = 1   
    evtwt_lhescale3 = 1   
    evtwt_lhescale4 = 1   
    evtwt_lhescale6 = 1   
    evtwt_lhescale8 = 1   
    evtwt_LHEScaleDownWeight  = 1
    evtwt_LHEScaleUpWeight = 1

    if (HEMwt == 0): 
      Jet_eta_phi_2D.Fill(entry.LeadJet_eta, entry.LeadJet_phi, evtwt)
      Lep_eta_phi_2D.Fill(entry.Lepton_eta, entry.Lepton_phi, evtwt)

    #Apply HEM
    #if (isMC(file) and year == "2018" and HEMwt == 0):
      #ran_proportion = random.uniform(0, 1)
      #print "ran_proportion", ran_proportion                                                                                               
      #if (ran_proportion > 0.65): HEMwt = 1
      #HEMwt = 0.35
    #evtwt  *= HEMwt

    if (isMC(file)):
      #random.seed(0)
      #nBTag  = 0
      for j in range(0, entry.nJet):
         if (Jet_pt(file, entry, j, "nominal") > 50 and abs(entry.Jet_eta[j]) < 2.4 and JetID(entry.Jet_jetId[j]) >= 1):  ##Changed
          jet4vec = TLorentzVector()
          jet4vec.SetPtEtaPhiM(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j], entry.Jet_phi[j], Jet_mass(file, entry, j, "nominal"))
          if (jet4vec.DeltaR(Lep4vec) > 0.5):
            
            #ran = random.uniform(0, 1)
            if (useDeepCSV): 
                if (btagDeepCSV(entry.Jet_btagDeepB[j], "Medium") >= 1): 
                  btagSF_evt_tree *= entry.Jet_btagSF_deepcsv[j]
                else:
                  btagSF_evt_tree *= (1 - min(1.0 , Jet_bTagEff(entry, j)*entry.Jet_btagSF_deepcsv[j]))/(1 - Jet_bTagEff(entry, j))
            else: 
                if (btagDeepFLV(entry.Jet_btagDeepFlavB[j], "Medium") >= 1): 
                  btagSF_evt_tree *= entry.Jet_btagSF_deepjet[j]
                  btagSF_evt_treeUp *= entry.Jet_btagSF_deepjet_up[j]
                  btagSF_evt_treeDown *= entry.Jet_btagSF_deepjet_down[j]
                else:
                  btagSF_evt_tree *= (1 - min(1.0 , Jet_bTagEff(entry, j)*entry.Jet_btagSF_deepjet[j]))/(1 - Jet_bTagEff(entry, j))
                  btagSF_evt_treeUp *= (1 - min(1.0 , Jet_bTagEff(entry, j)*entry.Jet_btagSF_deepjet_up[j]))/(1 - Jet_bTagEff(entry, j))
                  btagSF_evt_treeDown *= (1 - min(1.0 , Jet_bTagEff(entry, j)*entry.Jet_btagSF_deepjet_down[j]))/(1 - Jet_bTagEff(entry, j))
      #Define weights here 
      #muHPtSF = applyHighPtMuonSF (file, entry, Lep4vec.P(), Lep4vec.Eta())
      #evtwt *= muHPtSF
      topptweight=1
      topptweightUp=1
      topptweightDown=1
      LHEScaleUpWeight  = 1
      LHEScaleDownWeight  = 1
      if (doSys):
        lhe = []
        if (hasattr(entry, 'nLHEScaleWeight') and hasattr(entry, 'LHEScaleWeight')):
          for nScale in range(entry.nLHEScaleWeight):
            if not (nScale == 0 or nScale == 5 or nScale == 7):
              #print "LHEScale : ", nScale, ", weight ", entry.LHEScaleWeight[nScale]
              lhe.append(entry.LHEScaleWeight[nScale])
          lhe.sort()
          if (len(lhe)): 
            LHEScaleDownWeight = lhe[0]
            LHEScaleUpWeight = lhe[5]
          #print "LHEScaleDownWeight: ", LHEScaleDownWeight, ", LHEScaleUpWeight: ", LHEScaleUpWeight

      if (year  == "2016"):
          if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
            idsf=float((0.45*float(entry.Muon_ID_GH_SF[entry.Lepton_idx]))+(0.55*float(entry.Muon_ID_BCDEF_SF[entry.Lepton_idx])))
            isosf=np.add(np.multiply(0.45,entry.Muon_ISO_GH_SF[entry.Lepton_idx]), np.multiply(0.55,entry.Muon_ISO_BCDEF_SF[entry.Lepton_idx]))
            trigsf=np.add(np.multiply(0.45,entry.Muon_Trigger_GH_SF[entry.Lepton_idx]), np.multiply(0.55,entry.Muon_Trigger_BCDEF_SF[entry.Lepton_idx]))
            evtwt *= (float(idsf)*float(isosf)*float(trigsf))
          elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
            idisosf     = float(entry.Electron_IDISO_BCDEFGH_SF[entry.Lepton_idx])
            recosf      = float(entry.Electron_RECO_BCDEFGH_SF[entry.Lepton_idx])
            trigsf      = float(entry.Electron_Trigger_BCDEFGH_SF[entry.Lepton_idx])
            evtwt  *= (float(idisosf)*float(recosf)*float(trigsf))
          if (isTTJet(file) and applyTopPtRew): 
            #topptweight=float(entry.TopPt_Weight_Nominal)
            #Fit with M_lmj classification
            #topptweight= topPt_Reweighting(file, entry, 0.0615*2.11, 0.0005*0.514)
            #Fit with M_lj classification
            #topptweight= topPt_Reweighting(file, entry, 0.0615*1.333, 0.0005*0.0266)
            #Fit with M_lmj classification and only beta (v2)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.0976)
            #based on ST_v2 fit, only beta
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.99)
            #based on ST_v2 fit, norm and beta
            topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, 0.0005*0.877)
            topptweightDown= topptweight*topptweight
      elif (year == "2017"):
          if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
            idsf=float(entry.Muon_ID_BCDEF_SF[entry.Lepton_idx])
            isosf=float(entry.Muon_ISO_BCDEF_SF[entry.Lepton_idx])
            trigsf=float(entry.Muon_Trigger_BCDEF_SF[entry.Lepton_idx])
            evtwt *= (float(idsf)*float(isosf)*float(trigsf))
          elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
            idisosf     = float(entry.Electron_IDISO_BCDEF_SF[entry.Lepton_idx])
            recosf      = float(entry.Electron_RECO_BCDEF_SF[entry.Lepton_idx])
            #trigsf      = float(entry.Electron_Trigger_BCDEF_SF[entry.Lepton_idx])
            evtwt  *= (float(idisosf)*float(recosf))
          if (isTTJet(file) and applyTopPtRew): 
            #topptweight=float(entry.TopPt_Weight_Nominal)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*2.11, 0.0005*0.514)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*1.333, 0.0005*0.0266)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.0976)
            #based on ST_v2 fit
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.99)
            #based on ST_v2 fit, norm and beta
            topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, 0.0005*0.877)
            topptweightDown= topptweight*topptweight
      elif (year == "2018"):
          if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
            idsf=float(entry.Muon_ID_ABCD_SF[entry.Lepton_idx])
            isosf=float(entry.Muon_ISO_ABCD_SF[entry.Lepton_idx])
            #trigsf=float(entry.Muon_Trigger_ABCD_SF[entry.Lepton_idx])
            trigsfBefore=float(trig_SF(file, tree_muTrigBefore, entry.Lepton_pt, entry.Lepton_eta, "Muon"))
            trigsfAfter=float(trig_SF(file, tree_muTrigAfter, entry.Lepton_pt, entry.Lepton_eta, "Muon"))
            trigsf = float(((trigsfBefore*7.0) + (trigsfAfter*53.0))/60.0)
            #print "idsf ", idsf, " isosf ", isosf, " trigsf", trigsf
            evtwt *= (float(idsf)*float(isosf)*float(trigsf))
          elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
            idisosf     = float(entry.Electron_IDISO_ABCD_SF[entry.Lepton_idx])
            recosf      = float(entry.Electron_RECO_ABCD_SF[entry.Lepton_idx])
            #trigsf      = float(entry.Electron_Trigger_ABCD_SF[entry.Lepton_idx])
            evtwt  *= (float(idisosf)*float(recosf))
          if (isTTJet(file) and applyTopPtRew): 
            #topptweight=float(entry.TopPt_Weight_Nominal)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*2.11, 0.0005*0.514)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*1.333, 0.0005*0.0266)
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.0976)
            #based on ST_v2 fit
            #topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, -0.0005*0.99)
            #based on ST_v2 fit, norm and beta
            topptweight= topPt_Reweighting(file, entry, 0.0615*0.0, 0.0005*0.877)
            topptweightDown= topptweight*topptweight

      #evtwt_nominal = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*float(topptweight)*tm_WJets_Scaling(file, tree, (entry.ST_v2+HT_Central-entry.LeadJet_pt))
      evtwt_nominal = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*float(topptweight)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
      if (doSys):
        evtwt_PileupUp = evtwt*float(entry.genWeight)*float(entry.puWeightUp)*float(btagSF_evt_tree)*float(topptweight)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_PileupDown = evtwt*float(entry.genWeight)*float(entry.puWeightDown)*float(btagSF_evt_tree)*float(topptweight)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_BTagSFUp = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_treeUp)*float(topptweight)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_BTagSFDown = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_treeDown)*float(topptweight)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_topptweightUp = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*float(topptweightUp)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_topptweightDown = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*float(topptweightDown)*float(wpt_WJets_Scaling(file, entry))*st_TTJets_Scaling(file, entry)
        evtwt_LHEScaleDownWeight =  evtwt_nominal*LHEScaleDownWeight
        evtwt_LHEScaleUpWeight =  evtwt_nominal*LHEScaleUpWeight
        if (hasattr(entry, 'LHEScaleWeight')  and entry.nLHEScaleWeight): 
          evtwt_lhescale1 =  evtwt_nominal*entry.LHEScaleWeight[1]  
          evtwt_lhescale2 =  evtwt_nominal*entry.LHEScaleWeight[2]  
          evtwt_lhescale3 =  evtwt_nominal*entry.LHEScaleWeight[3]  
          evtwt_lhescale4 =  evtwt_nominal*entry.LHEScaleWeight[4]  
          evtwt_lhescale6 =  evtwt_nominal*entry.LHEScaleWeight[6]  
          evtwt_lhescale8 =  evtwt_nominal*entry.LHEScaleWeight[8]  


    #WJets Selection
    #if (MET_Filters(file, entry) >= 1 and entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1 and entry.LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4 and entry.DR_LepleadJet > 0.5 and nBTag == 0 and MET_pt(file, entry) > 60 and entry.ST_v2 > 500 and  HT_Central > 400 and MassT > 40 and abs(entry.DPHI_LepMet) < 0.5 and abs(entry.DPHI_LepleadJet) > 2.0 and (entry.nVetoMuons + entry.nVetoElectrons) == 1):
    cut = cutFlow(file, entry, "WJets", leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "nominal"), float(entry.LeadJet_pt))
    for itr in range (0, int(cut)):
      if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
        hmap["Mu_WJets_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          Mu_WJets_HT_vs_ST.Fill(HT_Central, entry.ST_v2, evtwt_nominal)
          fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt, "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,"topptweightUp", evtwt_topptweightUp)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,"topptweightDown", evtwt_topptweightDown)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt, "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt, "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("WJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)

      elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
        hmap["Ele_WJets_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          Ele_WJets_HT_vs_ST.Fill(HT_Central, entry.ST_v2, evtwt_nominal)
          fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightUp", evtwt_topptweightUp)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightDown", evtwt_topptweightDown)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("WJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)


    #Signal Pre Selection
    #BTag Eff Map after a preselection
    if (isMC(file)):
      if (MET_Filters(file, entry) >= 1 and entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1 and entry.LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4 and entry.DR_LepleadJet > 0.5 and MET_pt(file, entry, "nominal") > 60 and entry.ST_v2 > 500 and  HT_Central > 500 and MassT < 130):
        if (entry.nJets > 0):
          for j in range(0, entry.nJets):
            if (Jet_pt(file, entry, j, "nominal") > 30 and abs(entry.Jet_eta[j]) < 2.4 and JetID(entry.Jet_jetId[j]) >= 1):
              if (Jet_flavour(entry, j) == 0): effmap_pt_eta_b_all.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j])
              elif (Jet_flavour(entry, j) == 1): effmap_pt_eta_c_all.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j]) 
              elif (Jet_flavour(entry, j) == 2): effmap_pt_eta_l_all.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j]) 
          
              if (Jet_pt(file, entry, j, "nominal") > 30 and abs(entry.Jet_eta[j]) < 2.4 and JetID(entry.Jet_jetId[j]) >= 1 and btagDeepFLV(entry.Jet_btagDeepFlavB[j], "Medium") >= 1):
                if (Jet_flavour(entry, j) == 0): effmap_pt_eta_b_btagged.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j]) 
                elif (Jet_flavour(entry, j) == 1): effmap_pt_eta_c_btagged.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j]) 
                elif (Jet_flavour(entry, j) == 2): effmap_pt_eta_l_btagged.Fill(Jet_pt(file, entry, j, "nominal"), entry.Jet_eta[j]) 

    cut = cutFlow(file, entry, "PreSig", leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "nominal"), float(entry.LeadJet_pt))
    for itr in range (0, int(cut)):
      if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
        hmap["Mu_PreSig_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,"topptweightUp", evtwt_topptweightUp)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,"topptweightDown", evtwt_topptweightDown)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt, "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt, "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("PreSig", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)


        if (itr == 4 and doSkim):
          #Load tree containers
          Event_flag[0] = entry.Event_flag
          Lepton_pt[0] = entry.Lepton_pt
          Lepton_eta[0] = entry.Lepton_eta
          Lepton_phi[0] = entry.Lepton_phi
          Lepton_mass[0] = entry.Lepton_mass
          LeadJet_pt[0] = entry.LeadJet_pt
          LeadJet_eta[0] = entry.LeadJet_eta
          LeadJet_phi[0] = entry.LeadJet_mass
          LeadJet_btagDeepCSV[0] = entry.LeadJet_btag_DeepCSV
          LeadJet_btagDeepFLV[0] = entry.LeadJet_btag_DeepFLV
          LeadJet_btag[0] = leadjetBTag
          print leadjetBTag
          MET[0] = MET_pt(file, entry, "nominal")
          MT[0] = MassT
          HT[0] = HT_Central
          ST[0] = entry.ST
          ST_v2[0] = entry.ST_v2
          nVetoElectrons[0] = entry.nVetoElectrons
          nVetoMuons[0] = entry.nVetoMuons
          Mass[0] = entry.Mass
          Mass_v2[0] = entry.Mass_v2
          DR_LepClosestJet[0] = entry.DR_LepClosestJet
          DR_LepleadJet[0] = entry.DR_LepleadJet
          DPHI_LepleadJet[0] = entry.DPHI_LepleadJet
          DPHI_LepMet[0] = entry.DPHI_LepMet
          DPHI_JetMet[0] = entry.DPHI_JetMet
          Angle_MuZ_Jet[0] = entry.Angle_MuZ_Jet
          Angle_MuJet_Met[0] = entry.Angle_MuJet_Met
          bVsW_ratio[0] = entry.bVsW_ratio
          nJets[0] = entry.nJets
          #nCentralJets[0] = entry.nCentralJets
          nCentralJets[0] = nCentralJets_v2
          nFwdJets[0] = entry.nFwdJets
          nBTagMed_DeepCSV[0] = entry.nBTagMed_DeepCSV
          nBTagTig_DeepCSV[0] = entry.nBTagTig_DeepCSV
          nBTagMed_DeepFLV[0] = nBTag
          nBTagTig_DeepFLV[0] = entry.nBTagTig_DeepFLV
          FwdJetPt[0] = entry.FwdJetPt
          FwdJetEta[0] = entry.FwdJetEta
          FwdJetPhi[0] = entry.FwdJetPhi
          STScale_TT[0] = float(st_TTJets_Scaling(file, entry))
          WPTScale_WJ[0] = float(wpt_WJets_Scaling(file, entry))
          EvtWt[0] = evtwt_nominal

      elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
        hmap["Ele_PreSig_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightUp", evtwt_topptweightUp)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightDown", evtwt_topptweightDown)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("PreSig", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)


        if (itr == 4 and doSkim):
          #Load tree containers
          Event_flag[0] = entry.Event_flag
          Lepton_pt[0] = entry.Lepton_pt
          Lepton_eta[0] = entry.Lepton_eta
          Lepton_phi[0] = entry.Lepton_phi
          Lepton_mass[0] = entry.Lepton_mass
          LeadJet_pt[0] = entry.LeadJet_pt
          LeadJet_eta[0] = entry.LeadJet_eta
          LeadJet_phi[0] = entry.LeadJet_mass
          LeadJet_btagDeepCSV[0] = entry.LeadJet_btag_DeepCSV
          LeadJet_btagDeepFLV[0] = entry.LeadJet_btag_DeepFLV
          LeadJet_btag[0] = leadjetBTag
          MET[0] = MET_pt(file, entry, "nominal")
          MT[0] = MassT
          HT[0] = HT_Central
          ST[0] = entry.ST
          ST_v2[0] = entry.ST_v2
          nVetoElectrons[0] = entry.nVetoElectrons
          nVetoMuons[0] = entry.nVetoMuons
          Mass[0] = entry.Mass
          Mass_v2[0] = entry.Mass_v2
          DR_LepClosestJet[0] = entry.DR_LepClosestJet
          DR_LepleadJet[0] = entry.DR_LepleadJet
          DPHI_LepleadJet[0] = entry.DPHI_LepleadJet
          DPHI_LepMet[0] = entry.DPHI_LepMet
          DPHI_JetMet[0] = entry.DPHI_JetMet
          Angle_MuZ_Jet[0] = entry.Angle_MuZ_Jet
          Angle_MuJet_Met[0] = entry.Angle_MuJet_Met
          bVsW_ratio[0] = entry.bVsW_ratio
          nJets[0] = entry.nJets
          #nCentralJets[0] = entry.nCentralJets
          nCentralJets[0] = nCentralJets_v2
          nFwdJets[0] = entry.nFwdJets
          nBTagMed_DeepCSV[0] = entry.nBTagMed_DeepCSV
          nBTagTig_DeepCSV[0] = entry.nBTagTig_DeepCSV
          nBTagMed_DeepFLV[0] = nBTag
          nBTagTig_DeepFLV[0] = entry.nBTagTig_DeepFLV
          FwdJetPt[0] = entry.FwdJetPt
          FwdJetEta[0] = entry.FwdJetEta
          FwdJetPhi[0] = entry.FwdJetPhi
          STScale_TT[0] = float(st_TTJets_Scaling(file, entry))
          WPTScale_WJ[0] = float(wpt_WJets_Scaling(file, entry))
          EvtWt[0] = evtwt_nominal

      if ((entry.MuonTrigFlag >= 1 or entry.EleTrigFlag >= 1) and itr == 4 and doSkim): sktree.Fill()

    #Signal Selection
    #if (MET_Filters(file, entry) >= 1 and entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1 and entry.LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4 and entry.DR_LepleadJet > 0.5 and entry.DR_LepClosestJet > 1.5 and abs(entry.DPHI_LepMet) < 0.5 and entry.bVsW_ratio < 1.4 and abs(entry.DPHI_LepleadJet) > 2.0 and nBTag_30 >= 1 and MET_pt(file, entry) > 60 and entry.ST_v2 > 700 and  HT_Central > 500 and MassT < 130 and leadjetBTag >= 1 and (entry.nVetoMuons + entry.nVetoElectrons) == 1 and entry.FwdJetPt > 30 and abs(entry.FwdJetEta) > 2.4):

    cut = cutFlow(file, entry, "Signal", leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "nominal"), float(entry.LeadJet_pt))
    for itr in range (0, int(cut)):
      if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
        hmap["Mu_Signal_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("Signal", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)

      elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
        hmap["Ele_Signal_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("Signal", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)

    #TTJets Selection
    #if (MET_Filters(file, entry) >= 1 and entry.Lepton_pt > 40 and abs(entry.Lepton_eta) < 2.1 and entry.LeadJet_pt > 200 and abs(entry.LeadJet_eta) < 2.4 and entry.DR_LepClosestJet <  1.5 and nBTag >= 2 and MET_pt(file, entry) > 60 and entry.ST_v2 > 200 and leadjetBTag  >= 1):

    evtwt_TTJets_ST_Scaling = 1
    evtwt_def = 1
    evtwt_def2 = 1
    evtwt_alpha_up = 1
    evtwt_alpha_down = 1
    evtwt_beta_up = 1
    evtwt_beta_down = 1
    evtwt_beta_up2 = 1
    evtwt_beta_down2 = 1
    topptweight_def = 1
    topptweight_def2 = 1
    topptweight_alpha_up = 1
    topptweight_alpha_down = 1
    topptweight_beta_up = 1
    topptweight_beta_down = 1
    topptweight_beta_up2 = 1
    topptweight_beta_down2 = 1
    if (isMC(file)):
      if (isTTJet(file) and applyTopPtRew): 
        topptweight_def  = topPt_Reweighting(file, entry, 0.0615, 0.0005)
        topptweight_def2  = topPt_Reweighting(file, entry, 0.0, 0.0005)
        topptweight_alpha_up  = topPt_Reweighting(file, entry, 0.0615*2.0, 0.0005)
        topptweight_alpha_down  = topPt_Reweighting(file, entry, 0.0615/2.0, 0.0005)
        topptweight_beta_up  = topPt_Reweighting(file, entry, 0.0615, 0.0005*4.0)
        topptweight_beta_down  = topPt_Reweighting(file, entry, 0.0615, 0.0005*0.0)
        topptweight_beta_up2  = topPt_Reweighting(file, entry, 0.0, 0.0005*2.0)
        topptweight_beta_down2  = topPt_Reweighting(file, entry, 0.0, -0.0005*2.0)
      evtwt_def = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_def
      evtwt_def2 = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_def2
      evtwt_alpha_up = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_alpha_up
      evtwt_alpha_down = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_alpha_down
      evtwt_beta_up = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_beta_up
      evtwt_beta_down = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_beta_down
      evtwt_beta_up2 = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_beta_up2
      evtwt_beta_down2 = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*topptweight_beta_down2
      evtwt_TTJets_ST_Scaling = evtwt*float(entry.genWeight)*float(entry.puWeight)*float(btagSF_evt_tree)*st_TTJets_Scaling(file, entry)

    cut = cutFlow(file, entry, "TTJets", leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "nominal"), float(entry.LeadJet_pt))
    for itr in range (0, int(cut)):
      if (itr == 16):
        met4vec = ROOT.TLorentzVector()
        lep4vec = ROOT.TLorentzVector()
        jet4vec = ROOT.TLorentzVector()
        subjet4vec = ROOT.TLorentzVector()
        met4vec.SetPtEtaPhiM(MET_pt(file, entry, "nominal"), 0, MET_phi(file, entry, "nominal"), 0)
        lep4vec.SetPtEtaPhiM(entry.Lepton_pt, entry.Lepton_eta, entry.Lepton_phi, entry.Lepton_mass) 
        jet4vec.SetPtEtaPhiM(entry.LeadJet_pt, entry.LeadJet_eta, entry.LeadJet_phi, entry.LeadJet_mass)
        subjet4vec.SetPtEtaPhiM(Jet_pt(file, entry, 1, "nominal"), entry.Jet_eta[1], entry.Jet_phi[1], Jet_mass(file, entry, 1, "nominal"))
        Mlj_lead = (lep4vec + jet4vec).M() 
        Mlj_sublead = (lep4vec + subjet4vec).M() 
        top_score = 999
        Pt_lmj_lead = (met4vec + lep4vec + jet4vec).Pt() 
        Pt_lmj_sublead = (met4vec + lep4vec + subjet4vec).Pt() 
        if (Mlj_lead > 173.0 and Mlj_sublead > 173.0): top_score = 3
        elif (Mlj_lead < 173.0 and Mlj_sublead < 173.0): top_score = 2
        elif (Mlj_lead < 173.0 and Mlj_sublead > 173.0): top_score = 1
        elif (Mlj_lead > 173.0 and Mlj_sublead < 173.0): top_score = 0

      if (entry.Event_flag == 13 and entry.MuonTrigFlag >= 1 and (entry.nVetoMuons + entry.nVetoElectrons) == 1):
        hmap["Mu_TTJets_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightUp", evtwt_topptweightUp)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightDown", evtwt_topptweightDown)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)
            fillHisto("TTJets", "Mu", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "TTJets_ST_Scaling", evtwt_TTJets_ST_Scaling)

        if (itr == 16):
          if (top_score == 0): 
            hmap["Mu_TTJets_Pt_lmj_select"].Fill(Pt_lmj_sublead, evtwt_def)
            hmap["Mu_TTJets_Pt_lmj_select2"].Fill(Pt_lmj_sublead, evtwt_def2)
            hmap["Mu_TTJets_Pt_lmj_select_alpha_up"].Fill(Pt_lmj_sublead, evtwt_alpha_up)
            hmap["Mu_TTJets_Pt_lmj_select_alpha_down"].Fill(Pt_lmj_sublead, evtwt_alpha_down)
            hmap["Mu_TTJets_Pt_lmj_select_beta_up"].Fill(Pt_lmj_sublead, evtwt_beta_up)
            hmap["Mu_TTJets_Pt_lmj_select_beta_down"].Fill(Pt_lmj_sublead, evtwt_beta_down)
            hmap["Mu_TTJets_Pt_lmj_select_beta_up2"].Fill(Pt_lmj_sublead, evtwt_beta_up2)
            hmap["Mu_TTJets_Pt_lmj_select_beta_down2"].Fill(Pt_lmj_sublead, evtwt_beta_down2)
            Mu_TTJets_Mass_lead_sublead_bin0.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 1): 
            hmap["Mu_TTJets_Pt_lmj_select"].Fill(Pt_lmj_lead, evtwt_def)
            hmap["Mu_TTJets_Pt_lmj_select2"].Fill(Pt_lmj_lead, evtwt_def2)
            hmap["Mu_TTJets_Pt_lmj_select_alpha_up"].Fill(Pt_lmj_lead, evtwt_alpha_up)
            hmap["Mu_TTJets_Pt_lmj_select_alpha_down"].Fill(Pt_lmj_lead, evtwt_alpha_down)
            hmap["Mu_TTJets_Pt_lmj_select_beta_up"].Fill(Pt_lmj_lead, evtwt_beta_up)
            hmap["Mu_TTJets_Pt_lmj_select_beta_down"].Fill(Pt_lmj_lead, evtwt_beta_down)
            hmap["Mu_TTJets_Pt_lmj_select_beta_up2"].Fill(Pt_lmj_lead, evtwt_beta_up2)
            hmap["Mu_TTJets_Pt_lmj_select_beta_down2"].Fill(Pt_lmj_lead, evtwt_beta_down2)
            Mu_TTJets_Mass_lead_sublead_bin1.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 2):
            Mu_TTJets_Mass_lead_sublead_bin2.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 3):
            Mu_TTJets_Mass_lead_sublead_bin3.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          hmap["Mu_TTJets_Top_Score"].Fill(top_score, evtwt_def)
          hmap["Mu_TTJets_ST_v2_def2"].Fill(entry.ST_v2, evtwt_def2)
          hmap["Mu_TTJets_ST_v2_beta_up2"].Fill(entry.ST_v2, evtwt_beta_up2)
          hmap["Mu_TTJets_ST_v2_beta_down2"].Fill(entry.ST_v2, evtwt_beta_down2)

      elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1 and (entry.nVetoMuons + entry.nVetoElectrons) == 1):
        hmap["Ele_TTJets_Counter"].Fill(itr, evtwt_nominal)
        if (itr == 16):
          fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "nominal", evtwt_nominal)
          if (doSys):
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightUp", evtwt_topptweightUp)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "topptweightDown", evtwt_topptweightDown)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFUp", evtwt_BTagSFUp)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "BTagSFDown", evtwt_BTagSFDown)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupUp", evtwt_PileupUp)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "PileupDown", evtwt_PileupDown)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale1", evtwt_lhescale1)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale2", evtwt_lhescale2)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale3", evtwt_lhescale3)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale4", evtwt_lhescale4)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale6", evtwt_lhescale6)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScale8", evtwt_lhescale8)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleDownWeight", evtwt_LHEScaleDownWeight)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "LHEScaleUpWeight", evtwt_LHEScaleUpWeight)
            fillHisto("TTJets", "Ele", file, entry, itr, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "nominal"), MET_phi(file, entry, "nominal"), entry.LeadJet_pt,  "TTJets_ST_Scaling", evtwt_TTJets_ST_Scaling)


        if (itr == 16):
          if (top_score == 0): 
            hmap["Ele_TTJets_Pt_lmj_select"].Fill(Pt_lmj_sublead, evtwt_def)
            hmap["Ele_TTJets_Pt_lmj_select2"].Fill(Pt_lmj_sublead, evtwt_def2)
            hmap["Ele_TTJets_Pt_lmj_select_alpha_up"].Fill(Pt_lmj_sublead, evtwt_alpha_up)
            hmap["Ele_TTJets_Pt_lmj_select_alpha_down"].Fill(Pt_lmj_sublead, evtwt_alpha_down)
            hmap["Ele_TTJets_Pt_lmj_select_beta_up"].Fill(Pt_lmj_sublead, evtwt_beta_up)
            hmap["Ele_TTJets_Pt_lmj_select_beta_down"].Fill(Pt_lmj_sublead, evtwt_beta_down)
            hmap["Ele_TTJets_Pt_lmj_select_beta_up2"].Fill(Pt_lmj_sublead, evtwt_beta_up2)
            hmap["Ele_TTJets_Pt_lmj_select_beta_down2"].Fill(Pt_lmj_sublead, evtwt_beta_down2)
            Ele_TTJets_Mass_lead_sublead_bin0.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 1): 
            hmap["Ele_TTJets_Pt_lmj_select"].Fill(Pt_lmj_lead, evtwt_def)
            hmap["Ele_TTJets_Pt_lmj_select2"].Fill(Pt_lmj_lead, evtwt_def2)
            hmap["Ele_TTJets_Pt_lmj_select_alpha_up"].Fill(Pt_lmj_lead, evtwt_alpha_up)
            hmap["Ele_TTJets_Pt_lmj_select_alpha_down"].Fill(Pt_lmj_lead, evtwt_alpha_down)
            hmap["Ele_TTJets_Pt_lmj_select_beta_up"].Fill(Pt_lmj_lead, evtwt_beta_up)
            hmap["Ele_TTJets_Pt_lmj_select_beta_down"].Fill(Pt_lmj_lead, evtwt_beta_down)
            hmap["Ele_TTJets_Pt_lmj_select_beta_up2"].Fill(Pt_lmj_lead, evtwt_beta_up2)
            hmap["Ele_TTJets_Pt_lmj_select_beta_down2"].Fill(Pt_lmj_lead, evtwt_beta_down2)
            Ele_TTJets_Mass_lead_sublead_bin1.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 2):
            Ele_TTJets_Mass_lead_sublead_bin2.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          elif (top_score == 3):
            Ele_TTJets_Mass_lead_sublead_bin3.Fill(Mlj_lead, Mlj_sublead, evtwt_def)
          hmap["Ele_TTJets_Top_Score"].Fill(top_score, evtwt_def)
          hmap["Ele_TTJets_ST_v2_def2"].Fill(entry.ST_v2, evtwt_def2)
          hmap["Ele_TTJets_ST_v2_beta_up2"].Fill(entry.ST_v2, evtwt_beta_up2)
          hmap["Ele_TTJets_ST_v2_beta_down2"].Fill(entry.ST_v2, evtwt_beta_down2)


    #Systematics for JES and JER
    if (doSys):
      for cr in range(len(controlR)):
        cutjerUp = cutFlow(file, entry, str(controlR[cr]), leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "jerUp"), Jet_pt(file, entry, entry.LeadJet_idx, "jerUp"))
        cutjerDown = cutFlow(file, entry, str(controlR[cr]), leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "jerDown"), Jet_pt(file, entry, entry.LeadJet_idx, "jerDown"))
        cutjesUp = cutFlow(file, entry, str(controlR[cr]), leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "jesUp"), Jet_pt(file, entry, entry.LeadJet_idx, "jesUp"))
        cutjesDown = cutFlow(file, entry, str(controlR[cr]), leadjetBTag, HT_Central, nBTag, MassT, MET_pt(file, entry, "jesDown"), Jet_pt(file, entry, entry.LeadJet_idx, "jesDown"))
        CHL = ""
        if (float(entry.Event_flag) == 13 and float(entry.MuonTrigFlag) >= 1):
          CHL = "Mu"
        elif (entry.Event_flag == 11 and entry.EleTrigFlag >= 1):
          CHL = "Ele"
        if ('Mu' in CHL or 'Ele' in CHL):
          if (int(cutjerUp) == 17):   fillHisto(str(controlR[cr]), CHL, file, entry, int(cutjerUp)-1, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "jerUp"), MET_phi(file, entry, "jerUp"), Jet_pt(file, entry, entry.LeadJet_idx, "jerUp"),  "jerUp", evtwt_nominal)
          if (int(cutjerDown) == 17): fillHisto(str(controlR[cr]), CHL, file, entry, int(cutjerDown)-1, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "jerDown"), MET_phi(file, entry, "jerDown"), Jet_pt(file, entry, entry.LeadJet_idx, "jerDown"), "jerDown", evtwt_nominal)
          if (int(cutjesUp) == 17):   fillHisto(str(controlR[cr]), CHL, file, entry, int(cutjesUp)-1, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "jesUp"), MET_phi(file, entry, "jesUp"), Jet_pt(file, entry, entry.LeadJet_idx, "jesUp"), "jesUp", evtwt_nominal)
          if (int(cutjesDown) == 17): fillHisto(str(controlR[cr]), CHL, file, entry, int(cutjesDown)-1, nBTag, MassT, HT_Central, HEMwt, nCentralJets_v2, FwdJetPt_L, FwdJetPt_M, FwdJetPt_T, MET_pt(file, entry, "jesDown"), MET_phi(file, entry, "jesDown"), Jet_pt(file, entry, entry.LeadJet_idx, "jesDown"),"jesDown", evtwt_nominal)

  #canv.Modified()
  #canv.Update()
  outputfile.Write()
  outputfile.Close()
  print "Processing completed ... DONE"
