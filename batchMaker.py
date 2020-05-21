import sys, math, os, re
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('Path', '/uscms_data/d3/amodak/log/Condor/submit/Condorout',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "path to store"
)
options.parseArguments()

PATH = options.Path
print PATH

files_2016 = ["SingleMuon_2016B", "SingleElectron_2016B", "SingleMuon_2016C", "SingleElectron_2016C", "SingleMuon_2016D", "SingleElectron_2016D", "SingleMuon_2016E", "SingleElectron_2016E", "SingleMuon_2016F", "SingleElectron_2016F", "SingleMuon_2016G", "SingleElectron_2016G", "SingleMuon_2016H", "SingleElectron_2016H", "DYJetsToLL_M-50_HT-100to200_2016",
"DYJetsToLL_M-50_HT-200to400_2016",
"DYJetsToLL_M-50_HT-400to600_2016",
"DYJetsToLL_M-50_HT-600to800_2016", 
"DYJetsToLL_M-50_HT-800to1200_2016",
"DYJetsToLL_M-50_HT-1200to2500_2016",
"DYJetsToLL_M-50_HT-2500toInf_2016","WJetsToLNu_HT-100To200_2016", "WJetsToLNu_HT-200To400_2016", "WJetsToLNu_HT-400To600_2016", "WJetsToLNu_HT-600To800_2016", "WJetsToLNu_HT-800To1200_2016", "WJetsToLNu_HT-1200To2500_2016", "WJetsToLNu_HT-2500ToInf_2016", "ST_tW_top_2016", "ST_tW_antitop_2016", "ST_t-channel_top_2016", "ST_t-channel_antitop_2016", "TTJets_madgraphMLM_2016", "QCD_Pt_170to300_2016", "QCD_Pt_300to470_2016", "QCD_Pt_470to600_2016", "QCD_Pt_600to800_2016", "QCD_Pt_800to1000_2016", "QCD_Pt_1000to1400_2016", "QCD_Pt_1400to1800_2016", "QCD_Pt_1800to2400_2016", "QCD_Pt_2400to3200_2016", "QCD_Pt_3200toInf_2016", "TprimeBToBW_M-700_2016", "TprimeBToBW_M-800_2016", "TprimeBToBW_M-1200_2016", "TprimeBToBW_M-1300_2016", "TprimeBToBW_M-1400_2016", "TprimeBToBW_M-1500_2016", "TprimeBToBW_M-1600_2016", "TprimeBToBW_M-1700_2016"]

files_2017 = ["SingleMuon_2017B", "SingleElectron_2017B", "SingleMuon_2017C", "SingleElectron_2017C", "SingleMuon_2017D", "SingleElectron_2017D", "SingleMuon_2017E", "SingleElectron_2017E", "SingleMuon_2017F", "SingleElectron_2017F", "DYJetsToLL_M-50_HT-100to200_2017",
"DYJetsToLL_M-50_HT-200to400_2017",
"DYJetsToLL_M-50_HT-400to600_2017",
"DYJetsToLL_M-50_HT-600to800_2017", 
"DYJetsToLL_M-50_HT-800to1200_2017",
"DYJetsToLL_M-50_HT-1200to2500_2017",
"DYJetsToLL_M-50_HT-2500toInf_2017","WJetsToLNu_HT-100To200_2017", "WJetsToLNu_HT-200To400_2017", "WJetsToLNu_HT-400To600_2017", "WJetsToLNu_HT-600To800_2017", "WJetsToLNu_HT-800To1200_2017", "WJetsToLNu_HT-1200To2500_2017", "WJetsToLNu_HT-2500ToInf_2017", "ST_tW_top_2017", "ST_tW_antitop_2017", "ST_t-channel_top_2017", "ST_t-channel_antitop_2017", "TTJets_madgraphMLM_2017", 
"QCD_HT100to200_2017",
"QCD_HT200to300_2017",
"QCD_HT300to500_2017",
"QCD_HT500to700_2017",
"QCD_HT700to1000_2017",
"QCD_HT1000to1500_2017",
"QCD_HT1500to2000_2017",
"QCD_HT2000toInf_2017",
"TprimeBToBW_M-1200_2017", "TprimeBToBW_M-1300_2017", "TprimeBToBW_M-1400_2017", "TprimeBToBW_M-1500_2017", "TprimeBToBW_M-1600_2017", "TprimeBToBW_M-1700_2017", "TprimeBToBW_M-1800_2017", "TprimeBToBW_M-1900_2017", "TprimeBToBW_M-2000_2017", "TprimeBToBW_M-2100_2017", "TprimeBToBW_M-2200_2017", "TprimeBToBW_M-2300_2017", "TprimeBToBW_M-2400_2017"]

files_2018 = ["SingleMuon_2018A", "EGamma_2018A", "SingleMuon_2018B", "EGamma_2018B", "SingleMuon_2018C", "EGamma_2018C", "SingleMuon_2018D", "EGamma_2018D", "DYJetsToLL_M-50_HT-100to200_2018",
"DYJetsToLL_M-50_HT-200to400_2018",
"DYJetsToLL_M-50_HT-400to600_2018",
"DYJetsToLL_M-50_HT-600to800_2018", 
"DYJetsToLL_M-50_HT-800to1200_2018",
"DYJetsToLL_M-50_HT-1200to2500_2018",
"DYJetsToLL_M-50_HT-2500toInf_2018","WJetsToLNu_HT-100To200_2018", "WJetsToLNu_HT-200To400_2018", "WJetsToLNu_HT-400To600_2018", "WJetsToLNu_HT-600To800_2018", "WJetsToLNu_HT-800To1200_2018", "WJetsToLNu_HT-1200To2500_2018", "WJetsToLNu_HT-2500ToInf_2018", "ST_tW_top_2018", "ST_tW_antitop_2018", "ST_t-channel_top_2018", "ST_t-channel_antitop_2018", "TTJets_madgraphMLM_2018", 
"QCD_HT100to200_2018",
"QCD_HT200to300_2018",
"QCD_HT300to500_2018",
"QCD_HT500to700_2018",
"QCD_HT700to1000_2018",
"QCD_HT1000to1500_2018",
"QCD_HT1500to2000_2018",
"QCD_HT2000toInf_2018",
"TprimeBToBW_M-1200_2018", "TprimeBToBW_M-1300_2018", "TprimeBToBW_M-1400_2018", "TprimeBToBW_M-1500_2018", "TprimeBToBW_M-1600_2018", "TprimeBToBW_M-1700_2018", "TprimeBToBW_M-1800_2018", "TprimeBToBW_M-1900_2018", "TprimeBToBW_M-2000_2018", "TprimeBToBW_M-2100_2018", "TprimeBToBW_M-2200_2018", "TprimeBToBW_M-2300_2018", "TprimeBToBW_M-2400_2018"]


years = ["2016", "2017", "2018"]
for year in range(len(years)):
    print str(years[year])
    files =[]
    if (str(years[year]) == "2016"): files = "skimfiles_2016.txt"
    elif (str(years[year]) == "2017"): files = "skimfiles_2017.txt"
    elif (str(years[year]) == "2018"): files = "skimfiles_2018.txt"
    inputFile = open('batchDummy.py')
    jobid = 1
    with open(files, 'r') as f:
       for content in f:
          #content = f.readlines()
          inputFile = open('batchDummy.py')
          log = ""
          if ("WJetsToLNu_HT-100To200" in str(content)):
             outfile = "WJetsToLNu_HT-100To200_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-200To400" in str(content)):
             outfile = "WJetsToLNu_HT-200To400_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-400To600" in str(content)):
             outfile = "WJetsToLNu_HT-400To600_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-600To800" in str(content)):
             outfile = "WJetsToLNu_HT-600To800_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-800To1200" in str(content)):
             outfile = "WJetsToLNu_HT-800To1200_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-1200To2500" in str(content)):
             outfile = "WJetsToLNu_HT-1200To2500_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_HT-2500ToInf" in str(content)):
             outfile = "WJetsToLNu_HT-2500ToInf_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("DYJetsToLL_M-50_HT-200to400" in str(content)):
             outfile = "DYJetsToLL_M-50_HT-200to400_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("DYJetsToLL_M-50_HT-400to600" in str(content)):
             outfile = "DYJetsToLL_M-50_HT-400to600_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("DYJetsToLL_M-50_HT-600to800" in str(content)):
             outfile = "DYJetsToLL_M-50_HT-600to800_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("DYJetsToLL_M-50_HT-800to1200" in str(content)):
             outfile = "DYJetsToLL_M-50_HT-800to1200_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8" in str(content)):
             outfile = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTToHadronic_TuneCP5_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8" in str(content)):
             outfile = "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8" in str(content)):
             outfile = "TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8" in str(content)):
             outfile = "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8" in str(content)):
             outfile = "TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8" in str(content)):
             outfile = "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8" in str(content)):
             outfile = "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("ST_t-channel_antitop" in str(content)):
             outfile = "ST_t-channel_antitop_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("ST_t-channel_top" in str(content)):
             outfile = "ST_t-channel_top_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("ST_tW_antitop" in str(content)):
             outfile = "ST_tW_antitop_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("ST_tW_top" in str(content)):
             outfile = "ST_tW_top_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          elif ("WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8" in str(content)):
             outfile = "WJetsToLNu_TuneCP5_madgraphMLM-pythia8_"+str(years[year])+"_ext"+str(jobid)+".root"
             log  = outfile.split('.')[0]
             jobid += 1
          else:
             name = content.rsplit('/', 1)[1]
             log = name.split('.')[0]       

          outputFile = open('batch_'+str(log)+'.jdl', 'w')
          QUEUE = '1'
          EXE = 'run_'+str(log)+'.sh'
          for line in inputFile:
             line = line.replace('queue', QUEUE)
             line = line.replace('path', PATH)
             line = line.replace('exe', EXE)
             outputFile.writelines(line)
          inputFile.close()
          outputFile.close()
