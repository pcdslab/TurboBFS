# TurboBFS
TurboBFS is a highly scalable GPU-based set of top-down and bottom-up BFS algorithms in the language of linear algebra. These algorithms are applicable to unweighted, directed and undirected graphs represented by sparse adjacency matrices in the Compressed Sparse Column (CSC) format, and the transpose of the Coordinate (COO) format, which were equally applicable to direct and undirected graphs. 

The following BFS algorithms were implemented in TurboBFS: 

1. TurboBFS_COOC_TD:   top-down BFS algorithms for graphs represented by sparse adjacency matrices in the COOC format.
2. TurboBFS_CSC_BU:    bottom-up BFS algorithms for graphs represented by sparse adjacency matrices in the CSC format.
3. TurboBFS_CSC_TD:    top-down BFS algorithms for graphs represented by sparse adjacency matrices in the CSC format.
4. TurboBFS_CSC_TDBU:  top-down/bottom-up BFS algorithms for graphs represented by sparse adjacency matrices in the CSC format.
 
The TurboBFS algorithms were designed to process sparse adjacency matrices selected from the SuiteSparse Matrix Collection in the Matrix Market format. 

More details of the design, implementation and experimental results obtained with the TurboBFS algorithms are given in our paper cited below.
# Prerequisites
This software has been tested on the following dependences:
* CUDA 10.1
* gcc 8.4.0 
* Ubuntu 16.04.7

# Install
Series of instructions to install and compile the code:

1. Download the software, there is a self-contained folder for each BFS algorithm. The folder graphData contains some examples of Matrix Market files. 

git clone https://github.com/pcdslab/TurboBFS.git

2. Compile the code with the Makefile. Please set the library and include directories paths on the Makefile available for your particular machine. Example

$ cd TurboBFS/TurboBFS_CSC_TD

$ make
# Run
1. Update the path of the graphData folder in the corresponding .sh file
2. Add permision for the .sh file, example: chmod +x run_bfstdcsc.sh
3. run the experiments with: ./run_bfstdcsc.sh

# Publications

If you use this software please cite our paper:

Oswaldo Artiles, and Fahad Saeed, “TurboBFS: GPU Based Breadth-First Search (BFS) Algorithms in the Language of Linear Algebra”, 2021 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), May 17-21, 2021

# Acknowledgements
This research was supported by the National Science Foundations (NSF) under the Award Numbers CAREER OAC-1925960. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Science Foundation. We would also like to acknowledge the donation of a K-40c Tesla GPU and a TITAN Xp GPU from NVIDIA which was used for all the GPU-based experiments performed in our paper.

