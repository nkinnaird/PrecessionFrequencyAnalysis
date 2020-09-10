# 3-31-20: Script to hadd many root files together in stages, which improves the speed of hadding. The variable maxFiles can be adjusted accordingly, where the default is typically 20, but a value of 5 is better when hadding histograms with many folders/iterations.

outfile=$1
shift
files=$*

maxFiles=20
nFiles=$(echo $files | tr ' ' '\n' | wc -l)
subFiles=""

if [[ $nFiles -le $maxFiles ]]
then
    #echo "command is ..."
    hadd -f $outfile $files
else
    nSteps=$(python -c "from math import ceil; print int(ceil(float($nFiles)/$maxFiles))")
    echo "Merging $nFiles in $nSteps steps ..."

    for (( i=1; i<=$nSteps; i++)) 
    do
      lines=$(($i*$maxFiles))
      linesfromend=$maxFiles

      if [[ $lines -ge $nFiles ]]
	  then 
	  linesfromend=$(($maxFiles + $nFiles - $lines)) 
      fi
      
      partFiles=$(echo $files | tr ' ' '\n' | head -n $lines | tail -n $linesfromend)
      #partFiles=$(echo $files | tr ' ' '\n' | head -n $lines | tail -n $maxFiles)
      partOut=$outfile.part.$i
      
      #partFiles is a list of intermediate root files
      #echo "partFiles" $partFiles
      
      #partOut is the name of the intermediate hist
      #echo "partOut" $partOut

      hadd -f $partOut $(echo $partFiles | tr '\n' ' ')
      subFiles=$subFiles" "$partOut
      echo " "
    done

    numInterFiles=$(echo $subFiles | wc -w)
    echo "numInterFiles" $numInterFiles
    if [[ $numInterFiles -ge $maxFiles ]]
	then
	echo "one more merging step required..."
	nInterSteps=$(python -c "from math import ceil; print int(ceil(float($numInterFiles)/$maxFiles))")
	for (( j=1; j<=$nInterSteps; j++)) 
	  do
	  interLines=$(($j*$maxFiles))
	  interLinesfromend=$maxFiles
	  
	  if [[ $interLines -ge $numInterFiles ]]
	      then 
	      interLinesfromend=$(($maxFiles + $numInterFiles - $interLines)) 
	  fi
      
	  interPartFiles=$(echo $subFiles | tr ' ' '\n' | head -n $interLines | tail -n $interLinesfromend)
	  interPartOut=$outfile.subpart.$j
	  
	  hadd -f $interPartOut $(echo $interPartFiles | tr '\n' ' ')
	  interSubFiles=$interSubFiles" "$interPartOut

	  echo " "

	done

	echo "Merging intermediate files after second loop"
	hadd -f $outfile $interSubFiles
    else
	echo "Merging intermediate files"
	hadd -f $outfile $subFiles
      
    fi

fi

# line to get the directory within which this script is located
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# echo "Script dir is $SCRIPTDIR"

echo "Purging hadded file of duplicate cycles"
root -l -b -q $SCRIPTDIR/PurgeDuplicates.C+\(\"$outfile\"\)


# Remove the intermediate files
rm -f *.root.part.*
rm -f *.root.subpart.*
