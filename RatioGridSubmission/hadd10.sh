# 3-31-20: Similar to the default hadd.sh script, but this script does not add things together at the very end. It is used to take a set of N files and hadd them in groups of 10 (or whatever number), returning a total of N/10 files.
# This script is typically more useful when applied to Toy MC, allowing to gather a higher number of stats per root file.

outfile=$1
# shift
# files=$*
files=$2 # pass in a file list, not sequence of files

maxFiles=10
# nFiles=$(echo $files | tr ' ' '\n' | wc -l)
nFiles=$(less $files | wc -l)
subFiles=""

echo $nFiles

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

      partFiles=$(less $files | tr ' ' '\n' | head -n $lines | tail -n $linesfromend)      
      # partFiles=$(echo $files | tr ' ' '\n' | head -n $lines | tail -n $linesfromend)
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

    # numInterFiles=$(echo $subFiles | wc -w)
    # echo "numInterFiles" $numInterFiles
    # if [[ $numInterFiles -ge $maxFiles ]]
    # then
    # 	echo "one more merging step required..."
    # 	nInterSteps=$(python -c "from math import ceil; print int(ceil(float($numInterFiles)/$maxFiles))")
    # 	for (( j=1; j<=$nInterSteps; j++)) 
    # 	  do
    # 	  interLines=$(($j*$maxFiles))
    # 	  interLinesfromend=$maxFiles
    	  
    # 	  if [[ $interLines -ge $numInterFiles ]]
    # 	      then 
    # 	      interLinesfromend=$(($maxFiles + $numInterFiles - $interLines)) 
    # 	  fi
          
    # 	  interPartFiles=$(echo $subFiles | tr ' ' '\n' | head -n $interLines | tail -n $interLinesfromend)
    # 	  interPartOut=$outfile.subpart.$j
    	  
    # 	  hadd -f $interPartOut $(echo $interPartFiles | tr '\n' ' ')
    # 	  interSubFiles=$interSubFiles" "$interPartOut

    # 	  echo " "

    # 	done

    # 	echo "Merging intermediate files after second loop"
    # 	hadd -f $outfile $interSubFiles

    # else
    # 	echo "Merging intermediate files"
    # 	hadd -f $outfile $subFiles
      
    # fi

fi

# Remove the intermediate files
# rm -f *.root.part.*
# rm -f *.root.subpart.*
