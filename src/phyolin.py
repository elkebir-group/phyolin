from docplex.cp.model import CpoModel
import numpy as np
import pandas as pd
import os
import sys



def phyolin_model(B):
    
    ROWS, COLS = B.shape
    TOTAL = np.sum(B)
    
    mdl = CpoModel(name="Phyolin")
    
    c = mdl.integer_var_list(COLS, 1, COLS, "C")

    
    #t = mdl.integer_var_list(ROWS,1, taxa, "T")
    #m = mdl.integer_var_list(COLS, 1, mutclusters, "M")
    


    x = [[mdl.binary_var( name="B_hat" + str(r) + str(c)) for c in range(COLS)] for r in range(ROWS)]
    
    #the ordering of the columns must be different
    mdl.add(mdl.all_diff(c))

    #ensure 1's cannot change
    for i in range(ROWS):
        for j in range(COLS):
            v = B[i,j]
            if v > 0:
                mdl.add(x[i][j]==1)

    #constraint to ensure a linear phylogeny
    for k in range(ROWS):
        for i in range(COLS):
            for j in range(COLS):
                mdl.add(mdl.if_then(c[i] < c[j], x[k][j] <= x[k][i]))

    
    #objective function
    flips =  sum(sum(x[i][j] for i in range(ROWS)) for j in range(COLS)) - TOTAL
    mdl.add(mdl.minimize(flips))

    return mdl, x, c
        

def phyolin_cluster(B, taxa, mutations):
    ROWS, COLS = B.shape
    TOTAL = np.sum(B)
   
    mdl = CpoModel(name="Phyolin")
    
    c = mdl.integer_var_list(COLS, 1, COLS, "C")

    
    t = mdl.integer_var_list(ROWS,1, taxa, "T")
    m = mdl.integer_var_list(COLS, 1, mutations, "M")


    


    x = [[mdl.binary_var( name="B_hat" + str(r) + str(c)) for c in range(COLS)] for r in range(ROWS)]

    #the ordering of the columns must be different
    mdl.add(mdl.all_diff(c))

    #ensure 1's cannot change
    for i in range(ROWS):
        for j in range(COLS):
            v = B[i,j]
            if v > 0:
                mdl.add(x[i][j]==1)

    #constraint to ensure a linear phylogeny
    for k in range(ROWS):
        for i in range(COLS):
            for j in range(COLS):
                mdl.add(mdl.if_then(c[i] < c[j], x[k][j] <= x[k][i]))

    for i in range(ROWS):
        for j in range(i,ROWS):
            if(i != j):
                for k in range(COLS):
                    mdl.add(mdl.if_then(t[i]==t[j], x[i][k]==x[j][i]))


    
    for j in range(mutations):
        
        mdl.add(sum(m[k]==(j+1) for k in range(COLS))>= 1)
    
    for i in range(taxa):
        mdl.add(sum(t[k]==(i+1) for k in range(ROWS))>= 1)
    
    #objective function
    flips =  sum(sum(x[i][j] for i in range(ROWS)) for j in range(COLS)) - TOTAL
    mdl.add(mdl.minimize(flips))
    
    return mdl, x ,c 




def solve(mdl, time_limit=30):
    #print("\nSolving model....")
    msol = mdl.solve(TimeLimit = time_limit)
    
    obj = msol.get_objective_values()

    return msol, obj

def flips(msol,x, B):
    changes = []
    ROWS, COLS = B.shape
    for i in range(ROWS):
        for j in range(COLS):
            var = msol.get_value(x[i][j])
            if B[i,j] != var:
                change = "B" + str(i) + str(j)
                changes.append(change)
    return changes

def estFN(obj, B):
    row, col = B.shape
    Bsize = row*col
    fn = obj/row*col

    return fn

def convertFiletoMatrix(fname):
    #fname = "/home/leah/Documents/UIUC/Research/Phyolin/phyolin/data/test_data.csv"
    input_dat = pd.DataFrame(pd.read_csv(fname))
    #B = input_dat.sample(cells).to_numpy()
    B = input_dat.to_numpy()
    return(B)
    

def main(fname):
   
    B = convertFiletoMatrix(fname)


    print("Input Matrix:")
    # B = np.array([[1,1,0,0],[1,0,1,1],[0,0,1,0], [0,0,0,1],[0,0,0,1]])

    print("Generating Model....")
    mod,x ,c = phyolin_model(B)
    #mod, x, c = phyolin_cluster(B,t, m)
    print("Solving Model....")
    msol, obj = solve(mod)
    print("Solve Complete")
    print("Objective Value: " + str(obj[0]))
    var_flips = flips(msol,x,B)
    print("Flips:")
    for f in var_flips:
        print(f)
    print("estimated FN Rate: " + str(estFN(obj[0], B)))
#print(msol.print_solution())
#print(obj[0])
# print(msol)



if __name__ == "__main__":
    fname = sys.argv[1]
    #cells = sys.argv[2]
    t = 10
    m = 5
    main(fname)

# print(test_ex[:,[2,3,1,0]])