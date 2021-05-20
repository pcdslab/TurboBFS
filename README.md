# TurboBFS
A highly scalable GPU-based set of top-down and bottom-up BFS algorithms in the language of  linear algebra. These algorithms are applicable to unweighted, directed and undirected graphs represented by sparse adjacency matrices in the Compressed Sparse Column (CSC) format, and the transpose of the Coordinate (COO) format, which were equally applicable to direct and undirected graphs. We implemented the following BFS algorithms: 
1. TurboBFS_
# Prerequisites
This software has been tested on the following dependences:
* CUDA 10.1
* gcc 5.4.0 
* Ubuntu 16.04.5

# Install
Series of instructions to compile and run your code

1. Download the software

git clone https://github.com/pcdslab/TurboBFS.git

2. Compile the code with the Makefile. Please set the library and include directories paths on the Makefile available for your particular machine

$ cd GPU-SFFT

$ make
# Run

