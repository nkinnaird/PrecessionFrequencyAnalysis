# 3-31-20: Script to produce histograms on the grid from TTrees.
# Histograms are made with ClusterTreeToHistsOrdered.C, but the macro can be switched if desired.
# The input is a list of ROOT TTree clusters, in a .txt format where the input files are located somewhere on /pnfs/. These could be either production TTree clusters, or personal ones.
# Use as: . ratioConvertTreesToHistsOrdered.sh inputFileList.txt
# There are various grid related options which can be adjusted down below, such as the lifetime, memory usage, etc.
# There is various code throughout the script which deals with file names, in order to produce unique output file names for each individual job, and copy them properly back to /pnfs/.
# Running small test jobs on the grid is a good thing to do to make sure things are being done and used correctly.

#!bin/echo Please source instead

if [ $# -ne 1 ]; then
  echo "You must supply an input file list."
  return
fi

#Setup all the stuff we need here
source /grid/fermiapp/products/common/etc/setups.sh
setup jobsub_client
setup fife_utils

inputFile="$1"

#Make sure input file exists
if [ -f "$inputFile" ] && [ `cat "$inputFile" | wc -l` -gt 0 ]; then
  nFiles=`cat "$inputFile" | wc -l`
  echo "Processing $inputFile with $nFiles files..."
else
  echo "Input file list is required.  This was not found or has 0 entries"
  return
fi

#######################################

TARNAME="ClusterTreeToHistsTar.tar.gz"

read -p "Update tar file? (y/n) " yn
if [ "$yn" == "y" ]; then
  tar cfz $TARNAME ../HistMaking/ClusterTreeToHistsOrdered.C ../ratioMacroHeaders
fi

#Make sure tar file exists
if [ ! -f ${TARNAME} ]; then
  echo "Tar ball \"${TARNAME}\" is required (and must contain ClusterTreeToHistsOrdered.C).  This was not found."
  return
fi

#######################################


# Make output directory
read -p "Please enter output folder name: " folderName
pnfsOutDir=/pnfs/GM2/scratch/users/${USER}/RatioClusterHists/$folderName
if [ ! -d $pnfsOutDir ]; then
  mkdir -p $pnfsOutDir
  chmod -R g+w $pnfsOutDir
fi
echo "Output files will appear in $pnfsOutDir"

# Make sure that there's not already an output file that we'd overwrite
for file in `cat "$inputFile"`; do
  # Grab base file name
  # fileName=`basename $file`
  fileName=${file##*/}
  # Strip everything before second to last '_' to get unique ID
  uniqueID=${fileName#*_*_}
  # Remove .root from end
  uniqueID=${uniqueID%.*}
  if [ -f ${pnfsOutDir}/clusterHists_ordered_${uniqueID}.root ]; then 
    echo "ERROR: ${pnfsOutDir}/clusterHists_ordered_${uniqueID}.root already exists and would be overwritten."
    return
  fi
done

# Add xrootd business to the start of file names
rm -f xrootd_ClusterTreeFileList.txt && touch xrootd_ClusterTreeFileList.txt
for file in `cat "$inputFile"`; do
  # Strip /pnfs/ from file name and write into new file
  longFileName=${file#/pnfs/*}
  echo root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/${longFileName} >> xrootd_ClusterTreeFileList.txt
done

# Copy filelist, fcl & .so file to pnfs so we can get it from grid jobs
if [ -f ${pnfsOutDir}/xrootd_ClusterTreeFileList.txt ]; then
  rm -f ${pnfsOutDir}/xrootd_ClusterTreeFileList.txt 
fi
if [ -f ${pnfsOutDir}/${TARNAME} ]; then
  rm -f ${pnfsOutDir}/${TARNAME} 
fi
sleep 5
ifdh cp xrootd_ClusterTreeFileList.txt ${pnfsOutDir}/xrootd_ClusterTreeFileList.txt
ifdh cp ${TARNAME} ${pnfsOutDir}/${TARNAME}
echo "Copied xrootd_ClusterTreeFileList.txt & ${TARNAME} to $pnfsOutDir"

# File command - this is the one we'll want to run to get the right file
fileCmd='inFile=`sed "$(($PROCESS+1))q;d" xrootd_ClusterTreeFileList.txt`'

# Make script that we'll want to run on the grid
# Be sure to escape any variables that need to be evaluated at execution time
cat <<EOF > runGridJob.sh
#Setup stuff to copy files back and forth

# Some nodes have newer kernels
# if [[ $OS == "SL6" ]]; then
#   case `uname -r` in
#     3.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
#     4.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
#   esac
# fi

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
setup ifdhc

#Work out which sub-run we're processing - need to copy xrootdFileList here or can't sed it when it's on pnfs
if [ ! -z "\`ifdh ls ${pnfsOutDir}/xrootd_ClusterTreeFileList.txt\`" ]; then
  ifdh cp ${pnfsOutDir}/xrootd_ClusterTreeFileList.txt xrootd_ClusterTreeFileList.txt
else 
  echo "${pnfsOutDir}/xrootd_ClusterTreeFileList.txt does not exist"
fi
$fileCmd
fileName=\${inFile##*/}
uniqueID=\${fileName#*_*_}
uniqueID=\${uniqueID%.*}
echo "Processing file with unique ID: \$uniqueID"

#Copy root macro
if [ ! -z "\`ifdh ls ${pnfsOutDir}/${TARNAME}\`" ]; then
  ifdh cp ${pnfsOutDir}/${TARNAME} ./${TARNAME}
else 
  echo "${pnfsOutDir}/${TARNAME} does not exist"
fi

#Untar the tar-ball
tar xvzf ${TARNAME}

#Setup ROOT
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
setup root v6_12_04e -q e15:prof

#####################################################
# Do we need help from nova library_shim? Yes, we do!
#####################################################

# Just setup the shim even if we're on an onsite node

export PRODUCTS=\$PRODUCTS:/cvmfs/nova.opensciencegrid.org/externals
setup library_shim v03.03

#Found Nova shim or not?
if [ \$? -ne 0 ]; then 
   echo -e "\n!!!Library_shim NOT FOUND for site \$GLIDEIN_Site!!!\n"
else
   echo "Library_shim should be found for site \$GLIDEIN_Site"
fi

export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBRARY_SHIM_SL6_LIB_PATH
echo "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH"

#Get credentials ready to stream data with xrootd - for running locally, will cause issues on the grid
# kx509
# voms-proxy-init -debug -hours 72 -rfc -noregen -voms=fermilab:/fermilab/gm2/Role=Analysis

#Run root macro
root -l -b -q HistMaking/ClusterTreeToHistsOrdered.C+\(\"\${inFile}\"\)

if [ \$? -ne 0 ]; then 
   echo -e "\n Root macro unsuccessfully executed - exiting without attempting to copy back file. \n"
   exit 1
fi

#Copy output back
if [ -f clusterHists_ordered_\${uniqueID}.root ]; then 
  ifdh mv clusterHists_ordered_\${uniqueID}.root ${pnfsOutDir}/clusterHists_ordered_\${uniqueID}.root
else 
  echo "clusterHists_ordered_\${uniqueID}.root does not exist.  ls reads:"
  ls
  exit 1 # exit with a non-zero exit status code which I believe should be propogated such that the job is seen to fail by the grid and properly accounted for in summary reports
fi
EOF


if [ -f ${pnfsOutDir}/runGridJob.sh ]; then
  rm -f ${pnfsOutDir}/runGridJob.sh
fi
ifdh cp runGridJob.sh ${pnfsOutDir}/runGridJob.sh
while [ ! -f ${pnfsOutDir}/runGridJob.sh ]; do
  echo "${pnfsOutDir}/runGridJob.sh not found after ifdh cp.  Sleeping 5..."
  sleep 5
done


#jobsub_submit options

# USAGE_MODEL="DEDICATED,OPPORTUNISTIC,OFFSITE" # offsite will fail for some percentage of grid jobs - will need to sort out what to do in those cases in the future
USAGE_MODEL="DEDICATED,OPPORTUNISTIC"
# USAGE_MODEL="OFFSITE"
LIFETIME_OPT="16h"
MEMORYUSAGE="2GB" # this number might need to be changed depending on how many events/iterations are being ran - always test locally first and then update this - ~= 2200 MB + 600 * #iterations (when memory reduction is turned off)


#Submit grid job
jobsub_submit -N $nFiles -G gm2 --OS=SL6 --resource-provides=usage_model=${USAGE_MODEL} --expected-lifetime=${LIFETIME_OPT} --memory=${MEMORYUSAGE} --role=Analysis file://${pnfsOutDir}/runGridJob.sh
# jobsub_submit -N $nFiles -G gm2 --OS=SL6 --resource-provides=usage_model=${USAGE_MODEL} --expected-lifetime=${LIFETIME_OPT} --memory=${MEMORYUSAGE} --append_condor_requirements='(TARGET.HAS_CVMFS_gm2_opensciencegrid_org==true)' --append_condor_requirements='(TARGET.GLIDEIN_Site\ isnt\ \"Wisconsin\")' --role=Analysis file://${pnfsOutDir}/runGridJob.sh
rm -f runGridJob.sh
rm -f xrootd_ClusterTreeFileList.txt
