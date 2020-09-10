# Check that samweb and root are installed
command -v samweb >/dev/null 2>&1 || { echo >&2 "samweb not found. You want:"; echo "source /grid/fermiapp/products/common/etc/setups.sh"; echo "setup fife_utils"; return; }

# Get sam dataset number
if [ "$#" -ne 2 ]; then
  echo "SAM dataset and list of processed sub-runs required as arguments"
  return
fi
dataset=$1
filelist=$2
echo "Checking files in $filelist against dataset $dataset..."

# Get list of sub-runs in this dataset
samweb -e gm2 list-definition-files $dataset > SAMFileList.txt

# Get run & sub-runs for the SAMFileList 
rm -f SAMSubRuns.txt && touch SAMSubRuns.txt
for line in `cat SAMFileList.txt`; do
  run=${line##*_}
  run=${run%%.*}
  subrun=${line##*${run}.}
  subrun=${subrun%%.*}
  echo "$run $subrun" >> SAMSubRuns.txt
done
  
# Convert into format of "run*1000 + sub-run"
awk '{print $1*1000 + $2}' SAMSubRuns.txt > tmp.txt
sort tmp.txt > SAMSubRunsTmp.txt
mv SAMSubRunsTmp.txt SAMSubRuns.txt
rm -f tmp.txt

# Sort inputfile list
sort $filelist > tmp.txt
mv tmp.txt InputFileList.txt

# Output runs that are missing
echo "These run/sub-runs  were in the SAM dataset and not in $filelist:"
comm -23 SAMSubRuns.txt InputFileList.txt > MissingFilesTmp.txt
rm -f MissingFiles.txt && touch MissingFiles.txt
for line in `cat MissingFilesTmp.txt`; do 
  run=$(( line / 1000 ))
  subrun=$(( line-1000*run ))
  printf "%05d %05d\n" $run $subrun
  printf "%05d %05d\n" $run $subrun >> MissingFiles.txt
done

rm -f InputFileList.txt SAMFileList.txt SAMSubRuns.txt MissingFilesTmp.txt