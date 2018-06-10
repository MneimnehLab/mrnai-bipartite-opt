# mrnai-bipartite-opt
Software to predict the optimal interaction structure of multiple RNAs, under the bipartite interaction model

## Building the tool
Navigate to the `mrnai-bipartite-opt` directory, then simply run
```
make all
```
This will build mrnai-bipartite-opt as well as the RNAup software.

## Running the tool
To generate the optimal structure based on the complete bipartite graph approach (with even and odd RNAs specified) navigate to the `mrnai-bipartite-opt` directory and run 
```
./bip-run.sh <INPUT_FILE>
```
where `<INPUT_FILE>` is the path to your input file. Some sample input files are provided in the `mrnai-bipartite-opt/inputs' directory.

An input file consists of at least 3 lines:
Line 1: All RNA sequences, in 5' -> 3' direction, separated by `&`s.
Line 2: Any names that you can use for the RNAs,  separated by `&`s.
Line 3: Comma-separated indices of even RNAs, followed by `|`, followed by comma-separated indices of odd RNAs
Line 4 and onwards: Optional meta data.
Note that the indices of RNAs correspond to their ordering in line 1.


The output will be a set of windows.

As an example, try running
```
./bip-run.sh inputs/yeast_trunc
```
