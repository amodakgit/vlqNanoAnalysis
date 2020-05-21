#!/bin/bash                                                                                                                                                                                                                           
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
eval `scramv1 project CMSSW CMSSW_10_2_14`
cd CMSSW_10_2_14/src/
eval `scram runtime -sh`
echo "CMSSW: "$CMSSW_BASE

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
#let "sample=${1}+1"
echo "xrdcp root://cmseos.fnal.gov//store/user/amodak/toConndor/analyse_with_mva_batch.py"
xrdcp root://cmseos.fnal.gov//store/user/amodak/toConndor/analyse_with_mva_batch.py .
xrdcp root://cmseos.fnal.gov//store/user/amodak/toConndor/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root .
xrdcp root://cmseos.fnal.gov//store/user/amodak/toConndor/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root .
echo "run python analyzer"

python analyse_with_mva_batch.py --inputfile INFILE --datayear YEAR --outputfile OUTFILE



echo "xrdcp *.root root://cmseos.fnal.gov/OUTPUT/SAMPLE"


xrdcp -f OUTFILE root://cmseos.fnal.gov/OUTPUT/OUTFILE


rm *.root
rm *.py
cd ../../
rm -rf CMSSW_10_2_14

ls
echo "DONE!"
