# 3-31-20: Script to produce Toy MC histograms on the grid.
# This script will pass the script "submitRatioToyEvents.sh" to the grid where the contained code will run whichever Toy MC macro is desired.
# If new Toy MC code is written that one wishes to run on the grid, it needs to be included down below in the tarball.
# There are various grid related options which can be adjusted down below, such as the lifetime, memory usage, etc.
# Running small test jobs on the grid is a good thing to do to make sure things are being done and used correctly.

   
[ ! -r $PWD/submitRatioToyEvents.sh ] && echo "Unable to access submitRatioToyEvents.sh at $PWD !" && exit

TARNAME="RatioToyEventGeneratorTar.tar.gz"

read -p "Update tar file? (y/n) " yn
if [ "$yn" == "y" ]; then

  # tar cfz $TARNAME ../ToyMC/makeRatioToyHists.C ../ToyMC/makeToyHistFromTH2.C ../ToyMC/TimeEnergyPairsHist.root ../ratioMacroHeaders
  # tar cfz $TARNAME ../ToyMC/*.C  ../ToyMC/VW_MC/*.C ../ToyMC/ComparisonMC/OriginalComparison/*.C ../ToyMC/ComparisonMC/DifferentRandomizations/*.C ../ToyMC/ComparisonMC/WWoVWRand/*.C  ../ratioMacroHeaders
  # tar cfz $TARNAME ../ToyMC/*.C  ../ToyMC/VW_MC/*.C ../ToyMC/ComparisonMC/OriginalComparison/*.C ../ToyMC/ComparisonMC/DifferentRandomizations/*.C ../ToyMC/ComparisonMC/WWoVWRand/*.C ../ToyMC/ComparisonMC/DifferentMethods/*.C ../ToyMC/TimeEnergyPairsHist.root ../ratioMacroHeaders
  # tar cfz $TARNAME ../ToyMC/ComparisonMC/DifferentAnalyzers/makeAnalyzerHists.C ../ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/fornick_60h.root

  tar cfzv $TARNAME -C /gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/DifferentAnalyzers/ makeAnalyzerHists.C -C /gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/ . -C /gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/ . -C /gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/ .

fi


echo "Add folder tag information:"
read folderNameTag
SCRATCH_DIR=/pnfs/GM2/scratch/users/nbk228/RatioToyEvents/$folderNameTag

   if [ ! -d ${SCRATCH_DIR} ] ; then
      mkdir -p ${SCRATCH_DIR}
      chmod -R g+w ${SCRATCH_DIR}
   fi

ifdh cp $TARNAME ${SCRATCH_DIR}/$TARNAME

TARFILE=${SCRATCH_DIR}/${TARNAME}
TARFILE_OPT="-f ${TARFILE}"


MRB_PROJECT_VERSION="v9_21_06"

USAGE_MODEL="DEDICATED,OPPORTUNISTIC"
# USAGE_MODEL="DEDICATED,OPPORTUNISTIC,OFFSITE"

# LIFETIME_OPT="--expected-lifetime 16h"
LIFETIME_OPT="--expected-lifetime 24h"

MEMORYUSAGE="2GB" # this number actually needs to change depending on how many events/iterations are being ran - always test locally first and then update this


echo "How many jobs do you want to run?"
read NJOBS

    jobsub_submit \
    	     -G gm2 \
           -N $NJOBS \
           --OS=SL6 \
    	     ${TARFILE_OPT} \
           -e LOCAL_TARFILE_NAME=${TARNAME} \
           -e IFDH_CP_MAXRETRIES=5 \
           ${LIFETIME_OPT} \
           --memory=${MEMORYUSAGE} \
           --resource-provides=usage_model=${USAGE_MODEL} \
           file://$PWD/submitRatioToyEvents.sh ${SCRATCH_DIR} ${MRB_PROJECT_VERSION}
