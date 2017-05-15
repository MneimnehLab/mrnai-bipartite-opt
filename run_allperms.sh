INPUT_FILE=$1

if [ "$#" -lt 1 ]; then
	
	echo "Usage: "
	echo "./run.sh   relative/path/to/sequences_in"
	exit 1
fi

head -n 3 $INPUT_FILE | rnaup_weights/weights_rnaup
build/allperms -k 4 < $INPUT_FILE
