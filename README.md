# TurboBFS
A highly scalable GPU-based set of top-down and bottom-up BFS algorithms in the language of  linear algebra. These algorithms are applicable to unweighted, directed and undirected graphs represented by sparse adjacency matrices in the Compressed Sparse Column (CSC) format, and the transpose of the Coordinate (COO) format, which were equally applicable to direct and undirected graphs. We implemented the following BFS algorithms: 
1. TurboBFS_COOC_TD:   top-down BFS for graphs represented by sparse adjacency matrices in the COOC format.
2. TurboBFS_CSC_BU:    bottom-up BFS for graphs represented by sparse adjacency matrices in the CSC format.
3. TurboBFS_CSC_TD:    top-down BFS for graphs represented by sparse adjacency matrices in the CSC format.
4. TurboBFS_CSC_TDBU:  top-down/bootoom-up BFS for graphs represented by sparse adjacency matrices in the CSC format.
# Prerequisites
This software has been tested on the following dependences:
* CUDA 10.1
* gcc 8.4.0 
* Ubuntu 16.04.7

# Install
Series of instructions to compile and run your code

1. Download the software, there is a self-contained folder for each BFS algorirthm

git clone https://github.com/pcdslab/TurboBFS.git

2. Compile the code with the Makefile. Please set the library and include directories paths on the Makefile available for your particular machine

$ cd TurboBFS

$ make
# Run
1. Update the path of the graphData folder in the corresponding .sh file
2. Add permision for the .sh file, example: chmod +x run_bfstdcsc.sh
3. run the code with: ./run_bfstdcsc.sh

