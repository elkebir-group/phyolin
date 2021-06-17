"""
Helper script to convert simulated AML patient data to the input format for 
the deep learning method.
"""
import numpy as np 
import os



def convertMat(fname, fname_out, n, m):
    """
    Convert a simulated single-cell matrix to the input format for the deep learning comparison method 
    :param fname: the path to the simulated input data
    :param fname_out: the path to output file
    :param n: the number of rows in output matrix (padding is used if input rows < n)
    :param m: the number of columns in output matrix (padding is used if input columns < n)
    :return:
    """
    with open(fname) as fi:
        Binput = np.loadtxt(fi, dtype=np.int8, skiprows=1, delimiter=",")
        
    rows, cols = Binput.shape
    if rows < n:
        newrows =np.zeros((n - rows, cols), dtype = np.int8)
        Binput = np.vstack((Binput,newrows ))
    
    if cols < m:
        newcols = np.zeros((n,m-cols), dtype = np.int8)
        Binput = np.hstack((Binput, newcols))

    with open(fname_out, "w+") as predfile:
        np.savetxt(predfile,Binput, fmt="%01d")
 

input_dir = "../data/sim_data/input"
output_dir = "../data/sim_data/predict"
aml_files = os.listdir(input_dir)
n = 9300
m = 7
for a in aml_files:

  
        aml_in = os.path.join(input_dir, a)
        aml_out = os.path.join(output_dir, a)
        convertMat(aml_in, aml_out, n, m)

