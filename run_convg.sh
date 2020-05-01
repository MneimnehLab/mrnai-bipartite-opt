INPUT_FILE=$1

if [ "$#" -lt 1 ]; then
    
    echo "Usage: "
    echo "./run.sh   relative/path/to/sequences_in"
    echo "           [--wrap]"
    exit 1
fi

head -n 3 $INPUT_FILE | rnaup_weights/weights_rnaup -gu 1
build/converge -k 4 -gs 0 -gu 1  < $INPUT_FILE
