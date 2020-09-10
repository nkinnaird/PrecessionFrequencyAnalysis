listOfFiles=$1
output=list_of_runs_subruns.txt
# output=list_of_runs_subruns_fills.txt

#clear file
>| $output

for file in `cat $listOfFiles`; do

	# echo "Processing file: $file"

	grep "Proccessing FillID" $file | while read -r line; do
		fillID=${line##* }
		run=${fillID:0:5}
		subRun=${fillID:5:3}
		fill=${fillID:8:3}
		echo "$run" "$subRun" >> $output
		# echo "$run" "$subRun" "$fill" >> $output
	done

	less $output | sort -u -o $output

done



# ways to speed up maybe: 
#1. sort only once at the very end
#2. put in fill numbers as well and don't sort at all
