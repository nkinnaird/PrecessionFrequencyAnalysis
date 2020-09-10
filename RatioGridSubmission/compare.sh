inputFile=$1
# output=compare.txt

#clear file
# >| $output

	cat $inputFile | while read line; do
		run=${line:0:5}
		subRun=${line:6:3}
		compare=${run}_00${subRun}
		echo $compare
		# echo "$run_00$subRun"
	done


