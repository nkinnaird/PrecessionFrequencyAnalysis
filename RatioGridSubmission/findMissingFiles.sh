# 3-31-20: Bash script to look for missing files made with ClusterTreeToHistsOrdered.C - pass in the list of input files to the job and the output directory on /pnfs/.
# Can also be used to check for missing files made with LostMuonHistograms.C

if [ $# -ne 2 ]; then
  echo "You must supply two files (list of input files, output directory to search)."
  return
fi

inFile=$1
outDir=$2

for file in `cat $inFile`; do

  fileName=`basename $file`
  uniqueID=${fileName#*_*_}

  if [ ! -f $outDir/clusterHists_ordered_$uniqueID ]; then
  # if [ ! -f $outDir/lostMuons_hists_$uniqueID ]; then
    echo $file
  fi
done
