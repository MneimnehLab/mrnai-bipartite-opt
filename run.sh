INPUT_FILE=$1

if [ "$#" -lt 1 ]; then
    
    echo "Usage: "
    echo "./run.sh   relative/path/to/sequences_in"
    echo "           [--wrap]"
    exit 1
fi

head -n 3 $INPUT_FILE | rnaup_weights/weights_rnaup -gu 1
build/findone -k 4 -gs 0  < $INPUT_FILE
head -n 3 $INPUT_FILE | ./compile_weights.py > output/compiled_weights.out
