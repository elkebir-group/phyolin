# Phyolin - Identifying a linear perfect phylogeny in single-cell DNA sequencing data of tumors

The input to PhyDOSE is a binary matrix and a threshold false negative rate.  Under the assumption that the phylogeny is linear, Phyolin outputs an estimated false negative rate, an inferred linear phylogeny and classification of the tree topology based on the inputted threshold. 

![Overview of Phyolin](Figure1.png)

## Contents

  1. [Getting started](#start)
     * [Dependencies](#dep)
     * [Compilation](#comp)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
     * [Example](#example)

     

<a name="start"></a>
## Getting started

PhyDOSE is implemented in C++. 

| Folder    | DESCRIPTION                                                  |
| --------- | ------------------------------------------------------------ |
| `src`     | source code for PhyDOSE                                      |
| `data`    | example, simulated and real data for Phyolin                             
| `plots`   | output plots 


<a name="dep"></a>

### Dependencies   

PhyDOSE has the following dependencies:


* [CMake](http://www.cmake.org/) (>= 2.8)
* [Boost](http://www.boost.org) (>= 1.38)
* [CP Optimizer](https://www.ibm.com/analytics/data-science/prescriptive-analytics/cplex-optimizer) (>= 12.7)

<a name="comp"></a>
### Compilation

To compile Phyolin C++, execute the following commands from the root of the repository:

    $ mkdir build
    $ cd build
    $ cmake .. -DCPLEX=1
    $ make 


The compilation results in the following files in the `build` directory:

EXECUTABLE | DESCRIPTION
-----------|-------------
`phyolin`  | estimate the false negative rate under a null model of a linear topology

<a name="usage"></a>
## Usage Instructions

<a name="io"></a>
### I/O formats
The input to Phyolin is a .csv file that contains n +1 rows and m columns where n is the the number of single-cells and m is the number of mutations. The first row must contain the mutation names or ids. All other entries in the .csv file should be either 1 is mutation j is present in cell i an 0 otherwise. 

The output is a .csv file with the same number of rows and columns with the linear phylogeny inferred by Phyolin by flipping potential false negatives.




```
  Usage: 
    ./phyolin [--help|-h|-help] input output fn
  Where:
    --help|-h|-help
       Print a short help message
    input
       Input
    output
       Output
    fn
        a threshold false negative rate (float)
```
<a name="example"></a>
## Example
The following is an example of how to use Phyolin.

```
cd data
../build/phyolin test_data.csv test_data_out.csv 0.2

```
