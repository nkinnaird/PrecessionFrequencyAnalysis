# 3-31-20: Second script to produce Toy MC histograms on the grid, which is passed to grid nodes using the ratioToyEventsSetupAndSubmit.sh script.
# Adjust which Toy MC is ran by changing which macro is called down below.

######################################################
# If we don't have the g-2 CVMFS, then bail
######################################################	
if [ ! -d /cvmfs/gm2.opensciencegrid.org ]; then
  echo “CVMFS repo seems to not be present. Sleeping and then exiting.”
  sleep 100
  exit 1
fi

#############################
# environment variables
#############################

OUTPUT_DIR="$1"
RELEASE="$2"

####################################################
# create the environment job file
####################################################

# Some nodes have newer kernels
# if [[ $OS == "SL6" ]]; then
#   case `uname -r` in
#     3.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
#     4.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
#   esac
# fi

#/////////////////////////////////////////////////////////////////////////////////////

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup #TODO figure out if this is needed
setup ifdhc -z /cvmfs/fermilab.opensciencegrid.org/products/common/db

####################################################
# setup the g-2 software
####################################################
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
setup gm2 ${RELEASE} -q prof
setup ifdh_art v2_05_00 -q e15:prof:s63

#####################################################
# Do we need help from nova library_shim? Yes, we do!
#####################################################

# Just setup the shim even if we're on an onsite node

export PRODUCTS=$PRODUCTS:/cvmfs/nova.opensciencegrid.org/externals
setup library_shim v03.03

#Found Nova shim or not?
if [ $? -ne 0 ]; then 
   echo -e "\n!!!Library_shim NOT FOUND for site $GLIDEIN_Site!!!\n"
else
   echo "Library_shim found for site $GLIDEIN_Site"
fi


#####################################################
#####################################################

   LOCAL_TARFILE=${CONDOR_DIR_INPUT}/${LOCAL_TARFILE_NAME}
   tar -xzf ${LOCAL_TARFILE} >/dev/null

#/////////////////////////////////////////////////////////////////////////////////////

   # run the job

   # root -l -b ToyMC/makeRatioToyHists.C+
   # root -l -b ToyMC/ComparisonMC/WWoVWRand/makeRatioToyHistsWWoVWRand.C+
   # root -l -b ToyMC/ComparisonMC/WWoVWRand/makeRatioToyHistsWWoVWRandSimple.C+
   # root -l -b ToyMC/ComparisonMC/DifferentRandomizations/makeRatioToyHistsDiffRands.C+
   # root -l -b ToyMC/VW_MC/makeRatioToyHistsVW.C+
   # root -l -b ToyMC/VW_MC/makeRatioToyHistsVWOnly.C+

   # mv ratioToyHists.root ratioToyHists_${CLUSTER}.${PROCESS}.root
   # ifdh cp -D ratioToyHists_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}

   # root -l -b ToyMC/makeToyHistFromTH2.C+\(\"ToyMC/TimeEnergyPairsHist.root\"\)
   # mv histFromPairs.root histFromPairs_${CLUSTER}.${PROCESS}.root
   # ifdh cp -D histFromPairs_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}


   # root -l -b ToyMC/ComparisonMC/DifferentMethods/makeToyHistsAlexModified.C+\(\"ToyMC/TimeEnergyPairsHist.root\"\)
   # mv toyHistsMethods.root toyHistsMethods_${CLUSTER}.${PROCESS}.root
   # ifdh cp -D toyHistsMethods_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}

   # root -l -b ToyMC/ComparisonMC/DifferentMethods/makeToyHistsNoRandomization.C+\(\"ToyMC/TimeEnergyPairsHist.root\"\)
   # mv toyHistsMethodsNoRand.root toyHistsMethodsNoRand_${CLUSTER}.${PROCESS}.root
   # ifdh cp -D toyHistsMethodsNoRand_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}

   # root -l -b ToyMC/ComparisonMC/DifferentAnalyzers/makeAnalyzerHists.C+\(\"ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/fornick_60h.root\"\)
   # mv analyzerHists.root analyzerHists_${CLUSTER}.${PROCESS}.root
   # ifdh cp -D analyzerHists_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}

   root -l -b makeAnalyzerHists.C+\(1\)
   mv analyzerHists.root analyzerHists_${CLUSTER}.${PROCESS}.root
   ifdh cp -D analyzerHists_${CLUSTER}.${PROCESS}.root ${OUTPUT_DIR}

#/////////////////////////////////////////////////////////////////////////////////////
  
   # for copying log/env files - had a problem with this last time I tried to uncomment it

   # ifdh cp -D job_output_${CLUSTER}.${PROCESS}.log ${OUTPUT_DIR}
   # ifdh cp -D env_${CLUSTER}.${PROCESS}.log ${OUTPUT_DIR}

