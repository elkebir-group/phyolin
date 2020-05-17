import numpy as np 
import os

def run_ms(filename, nCell, nMut):
    msAddress = "/scratch/data/leah/external/ms"
    cmd = f"{msAddress} {nCell} 1 -s {nMut} | tail -n {nCell} > {filename}"
    os.system(cmd)
    with open(filename, 'r') as f:
        l = [line for line in f]
    l1 = [s.strip('\n') for s in l]
    l2 = np.array([[int(s) for s in q] for q in l1])  # Original matrix
    return l2

def convertMat(fname):
    with open(fname) as fi:
        foo= np.loadtxt(fi, dtype=np.int8, skiprows=2)
        # skiprows=2)
    
    return(foo)
                

def addNoise(B, alpha=0.1, beta=0.2, delta=0.1):
    rows, cols = B.shape
    fp = 0
    fn = 0
    for i in range(rows):
        for j in range(cols):
            if B[i,j]==0 and np.random.random_sample() < alpha:
                B[i,j] = 1
                fp = fp + 1
            
            if B[i,j] == 1 and np.random.random_sample() < beta:
                B[i,j] = 0
                fn = fn + 1
    
    for i in range(rows):
        if np.random.random_sample() < delta:
            k = np.random.randint(0, rows-1, size =1)
            for j in range(cols):
                if B[k,j] ==1:
                    B[i,j] =1
    
    print(fp)
    print(fn)
    
    return B

def saveTrue(B, base_path, sim_id, label):
    rows, cols = B.shape
    fname = "sim" + str(sim_id) + "_m" + str(rows) + "_n" + str(cols) + "_" + label
    

    fname_csv = base_path + "/true/" + fname + ".csv"
    
    Bcsv = np.vstack((np.array(range(0, cols)), B))
    np.savetxt(fname_csv, Bcsv,fmt="%01d", delimiter="," )


def saveB(B, base_path, sim_id, label):
    rows, cols = B.shape
    fname = "sim" + str(sim_id) + "_m" + str(rows) + "_n" + str(cols) + "_" + label
    
    fname_txt = base_path + "/predict/" + fname + ".B"
    fname_csv = base_path + "/input/" + fname + ".csv"
    np.savetxt(fname_txt, B, fmt="%01d")
    
    Bcsv = np.vstack((np.array(range(0, cols)), B))
    np.savetxt(fname_csv, Bcsv,fmt="%01d", delimiter="," )
#take from https://github.com/alreadydone/PhyloM/blob/master/Branching_Inference/util.py
def make_constrained_linear_input(n, m, w):
  """
  :param n:
  :param m:
  :param w:  Number of ones, i.e., weight.
  :return:
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

#take from https://github.com/alreadydone/PhyloM/blob/master/Branching_Inference/util.py
def is_linear(x):
  colOrder = np.argsort(-np.sum(x, axis=0), )
  x = x[:, colOrder]
  y = x[np.lexsort(np.rot90(x))]
  for j in range(y.shape[1]-1):
    u = int(y.shape[0] - np.sum(y[:,j]))
    remain = np.count_nonzero(y[:u, (j+1):])
    # print(j, u, remain)
    # print(y[(j+1):, :u])
    if remain!= 0:
      return False
  return True



def genSims(base_path, n=100, m=100, num_sims = 25, 
                alpha= 0.002, beta=0.2, delta=0.1):
    linear_count = 0
    branch_count = 0
    labels = []
    sim_id = 0
    while(linear_count < num_sims and branch_count < num_sims):
        B =run_ms("temp.txt", m, n)
        if(is_linear(B)):
            linear_count= linear_count + 1
            Bnoisy = addNoise(B,  alpha= alpha, beta = beta, delta = delta)
            sim_id = sim_id + 1
            label = "linear"
            saveB(Bnoisy,base_path, sim_id, label )
            saveTrue(B, base_path, sim_id, label)
            labels.append(label)
        else:
            branch_count = branch_count + 1
            sim_id = sim_id + 1
            Bnoisy = addNoise(B,  alpha= alpha, beta = beta, delta = delta)
            label= "branched"
            labels.append(label)
            saveB(Bnoisy,base_path, sim_id, label )
            saveTrue(B, base_path, sim_id, label)
            weight = np.sum(B)

            B_lin = make_constrained_linear_input(n, m, weight)
            
            Bnoisy = addNoise(B_lin,  alpha= alpha, beta = beta, delta = delta)
            sim_id = sim_id + 1
            linear_count = linear_count + 1
            label = "linear"
            labels.append(label)
            saveTrue(B_lin, base_path, sim_id, label)
            saveB(Bnoisy, base_path, sim_id, label )



#genSims("../data/sim2", num_sims=25, delta=0)
#genSims("../data/sim3", num_sims=25, delta=.1)
genSims("../data/sim4", num_sims=25, beta = 0.4, delta = 0.2)
genSims("../data/sim5", num_sims=25, beta = 0.4, delta = 0.0)
genSims("../data/sim6", num_sims=25, beta = 0.4, delta = 0.1)
genSims("../data/sim7", num_sims=25, beta = 0.4, delta = 0.2)

#print(is_linear(foo))
#print(is_linear(make_constrained_linear_input(100,100, np.sum(foo))))
#fname = "../data/sim2/test.B"
#myB = convertMat(fname)

#noisyB = addNoise(myB)
#saveB(noisyB, "../data/sim2/noisytest.B")
# print(make_constrained_linear_input(100,100, 20))


