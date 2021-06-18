#!/usr/bin/env python3
"""Generate simulated single-cell data with errors (false negatives, false positives and doublets) 
that represents both linear and branched phylogenies"""

import numpy as np 
import os
import sys
import argparse

def run_ms(filename, nCell, nMut):
    """
    Call the ms program to generate error free single-cell data that represent a
    perfect phylogeny. 
    :param filename: the path where the ms generated error-free single-cell data should be written
    :param nCell: the number of cells
    :param nMut: the number of mutations
    :return: an numpy array (shape: nCell x nMuts) of error free binary single-cell data
    """
    msAddress = "/scratch/data/leah/external/ms"
    cmd = f"{msAddress} {nCell} 1 -s {nMut} | tail -n {nCell} > {filename}"
    os.system(cmd)
    with open(filename, 'r') as f:
        l = [line for line in f]
    l1 = [s.strip('\n') for s in l]
    l2 = np.array([[int(s) for s in q] for q in l1])  
    return l2

                

def addNoise(B, alpha=0.0, beta=0.0, delta=0.0, missing=0):
    """
    Introduce techinical errors to error free single-cell data
    :param B: a numpy array of error free single-cell data
    :param alpha: the false positive rate
    :param beta: the false negative rate
    :param delta: the doublet rate
    :param missing: the missing data rate
    :return: an numpy array (shape: nCell x nMuts) of single-cell data with errors, 
    the count of simulated false positives, the count of simulated false negatives and 
    the count of missing entries.
    """
    rows, cols = B.shape
    Bout = B.copy()
    fp = 0
    fn = 0
    missing_count = 0
    for i in range(rows):
        for j in range(cols):
            if B[i,j]==0 and np.random.random_sample() < alpha:
                Bout[i,j] = 1
                fp = fp + 1
            
            if B[i,j] == 1 and np.random.random_sample() < beta:
                Bout[i,j] = 0
                fn = fn + 1
            
            if np.random.random_sample() < missing:
                Bout[i,j] == -1
                missing_count = missing_count + 1
    
    for i in range(rows):
        if np.random.random_sample() < delta:
            k = np.random.randint(0, rows-1, size =1)
            for j in range(cols):
                if B[k,j] ==1:
                    Bout[i,j] =1
    
    
    return Bout, fp, fn, missing_count

def saveTrue(Bin, fname):
    """
    Save the simulated error-free matrix single-cell data 
    :param Bin: a numpy array of the error free simulated single-cell data
    :param fname: the path to output file
    :return: 
    """
    rows, cols = Bin.shape
    Bcsv = np.vstack((np.array(range(0, cols)), Bin))
    np.savetxt(fname, Bcsv,fmt="%01d", delimiter="," )


def saveB(B, fname):
    """
    Save the simulated single-cell data in a Phyolin ready format
    :param B: a numpy array of the simulated single-cell data
    :param fname: the path to output file
    :return: 
    """
    np.savetxt(fname, B, fmt="%01d")
    


def savePhyInput(B, fname):
    """
    Save the simulated single-cell data in a Phyolin ready format
    :param B: a numpy array of the simulated single-cell data
    :param fname: the path to output file
    :return: 
    """
    rows, cols = B.shape
    print(f"row:{rows} col:{cols}")
    
    np.savetxt(fname, B,fmt="%01d", delimiter="," )


def make_constrained_linear_input(n, m, w):
  """
  Convert a branched phylogeny to a linear phylogeny
  :param n: the number of cells
  :param m: the number of mutations
  :param w:  Number of ones, i.e., weight.
  :return: a numpy array of single-cell data that represents linear phylogeny
  """
  x = np.zeros((n, m))
  u = 0
  for j in range(m):
    if n - w > u:
      u = n - w
    uu = int(n - w / (m - j))+1
    u = np.random.randint(u, uu)
    w -= n-u
    x[u:, j] = 1
  return x

def is_linear(x):
  """
  Check if the simulated error-free single-cell matrix represents a linear phylogeny
  :param x: a numpy array of the error-free single-cell data
  :return: a boolean indicating if the input matrix represents a linear phylogeny
  """
  colOrder = np.argsort(-np.sum(x, axis=0), )
  x = x[:, colOrder]
  y = x[np.lexsort(np.rot90(x))]
  for j in range(y.shape[1]-1):
    u = int(y.shape[0] - np.sum(y[:,j]))
    remain = np.count_nonzero(y[:u, (j+1):])
    if remain!= 0:
      return False
  return True




def main(args):

  
    if args.type=="linear":
        phy_type = 0  
    else:
        phy_type = 1

    n = args.cells
    m = args.mutations
    if m < 0 or n < 0:
        raise Exception("Desired cells and loci must be greater than 0")


    true_file = args.true
    alpha = args.alpha
    delta = args.delta
    beta = args.beta 
    outfile = args.outputfile
    countfile = args.counts
    predfile = args.predict
    missing = args.missing

    print(f"missing:{missing}")
    
    while(True):
        B =run_ms(args.msfile, nCell=n, nMut=m)
        print(B.shape)
        is_lin = is_linear(B)
        if not is_lin and phy_type ==0:
            B = make_constrained_linear_input(n, m,np.sum(B))
            break

         
        elif is_lin and phy_type==0:
            break
        elif not is_lin and phy_type ==1:
            break

        else:
            continue


      
    saveTrue(B, true_file)
    Bnoisy, fp, fn, missing_count = addNoise(B,  alpha= alpha, beta = beta, delta = delta, missing=missing)
    saveB(Bnoisy, predfile)      
    savePhyInput(Bnoisy, outfile )
    with open(countfile, 'w') as output:
            output.write("label\tcells\tloci\tfp\tfn\tmissing\n")
            
            output.write(f"{args.type}\t{n}\t{m}\t{fp}\t{fn}\t{missing_count}\n")

            


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--mutations", type=int, default=100, help="number of mutations")
    parser.add_argument("-n", "--cells", type=int, default=100, help="number of cells")
    parser.add_argument("--beta", type=float, default=0.0, help="Allelic dropout (ADO) rate [0.05]")
    parser.add_argument("--delta", type=float, default=0, help="doublet rate [0.1]")
    parser.add_argument("--alpha", type=float, default = 0, help="false positive rate [0]")
    parser.add_argument("--missing", type=float, default = 0, help="fraction of data this missing [0]")
    parser.add_argument("--type", type=str, default = "linear", help="type, one of linear or branched")
    parser.add_argument("--predict", type=str, help="type, one of linear or branched")
    parser.add_argument("--true", type=str, help="true phylogeny filename")
    parser.add_argument("--counts", type=str, help="record of fp and fn counts filename")
    parser.add_argument("--msfile", type=str, help="ms filename")
    parser.add_argument("-o", "--outputfile", type=str, help="output file name")
    

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)




