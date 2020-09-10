inFile=$1

cat $inFile | grep "IFFileTransfer: entered for uri:" | awk -F "usr" '{print $2}' | awk -v prefix="/pnfs" '{print prefix $0}'
