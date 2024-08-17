# -*- coding: utf-8 -*-

"""
Resampling maps program
"""

from random import random as rn
import sys
import numpy as np
from scipy.stats import norm
#from scipy.linalg import lu
from copy import deepcopy
import pandas as pd
from numpy import linalg as la

import numpy as np #,numpy.linalg

from datetime import datetime


rep = 100


#%% Functions etc
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

def pd(B):
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False


def nearpd(A):
    B = (A + A.T) / 2
    _, s, V = la.svd(B)
    H = np.dot(V.T, np.dot(np.diag(s), V))
    A2 = (B + H) / 2
    A3 = (A2 + A2.T) / 2

    if pd(A3):
        return A3
    
    spacing = np.spacing(la.norm(A))
    
    I = np.eye(A.shape[0])
    k = 1
    while not pd(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1
    
    return A3


#%% END function etc

#%% File names

# Input file ------------------------------------------------------------------
file_c = r"out_corr.csv"         # Correlation matrix file
file_p = r"out_par.csv"          # Parameters and std dev file
file_s = r"out_std_I_changed_it.csv"
file_s = r"out_std_from_correlation.csv"
path     = r"D:\work" + "\\"


# Output file -----------------------------------------------------------------
path_out = r"D:\work" + "\\"
file_out = r"output_samples_STD_FROM_CORRELATION.out"

#%% END of File names

par = pd.read_csv(path + file_p).values.tolist()
par = np.array(par[0])

std = pd.read_csv(path + file_s).values.tolist()
std = np.array(std[0])

corr = pd.read_csv(path + file_c).values.tolist()
corr = np.array(corr)


print("len  par", len(par))
print("len  std", len(std))
print("len corr", len(corr))


print("starting Nearest Corr")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

par =  np.array(par)
std =  np.array(std)
corr = np.array(corr)
M, N = corr.shape            # lines M, columns N ... vector = line
diag = np.zeros((M, N))      # matrix with diagonals as std errors
for i in range(M):           # loop to make corr matrix symmetric and create diagonal matrix with std error
    diag[i][i] = std[i]
    for ii in range(N):
        corr[i][ii] = corr[ii][i]

# Force positive definite
corr = deepcopy( nearpd(corr) )


#corr = deepcopy(nearPSD(corr ))


print("starting Multiplications")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

### Declaring variables, matrices, etc ----------------------------
n_par = len(par)             # the number of parameters
#cov = np.zeros((M, N))       # covariance matrix
Z = np.zeros(n_par)
tmp_samp = np.zeros(n_par)   # temporary sampled parameters

samp_par = [0 for i in range(rep)] # create list of np.zeros with M,N dimension
### Endo of Declaration -------------------------------------------
cov = np.matmul( np.matmul(diag, corr), diag )
#cov = nearpd(cov)

print("Starting Cholesky")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

L = la.cholesky(cov)


for i in range(rep):
    if 1 % 10 == 0:
        print("repetition:", i + 1)
    for iii in range(n_par):
        Z[iii] = norm.ppf(rn())  # create vector with random percentiles of std normal dist
    tmp_samp = np.matmul(L, Z) + par
    samp_par[i] = deepcopy(tmp_samp)


df = pd.DataFrame(samp_par)
df.to_csv(path_out + file_out)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)



sys.exit()
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END
# END END END END END END END END END END END END END END END END END END END END END END END END END END



def cho(par = par, std = std, corr = corr, rep = rep, M = M, N = N):
    
    print("starting cho")
    ### Declaring variables, matrices, etc ----------------------------
    n_par = len(par)             # the number of parameters
    #cov = np.zeros((M, N))       # covariance matrix
    Z = np.zeros(n_par)
    tmp_samp = np.zeros(n_par)   # temporary sampled parameters
    
    samp_par = [0 for i in range(rep)] # create list of np.zeros with M,N dimension
    ### Endo of Declaration -------------------------------------------
    cov = np.matmul( np.matmul(diag, corr), diag )
    L = la.cholesky(cov)
    
    
    for i in range(rep):
        print("repetition:", i)
        for iii in range(n_par):
            Z[iii] = norm.ppf(rn())  # create vector with random percentiles of std normal dist
        tmp_samp = np.matmul(L, Z) + par
        samp_par[i] = deepcopy(tmp_samp)
    
    return samp_par

aaa = cho()








# Function
def mapp(x,y):
    return x ** 0.5 + y ** 0.6


# Creating image
data = [ [0] * 11 for _ in range(11)]
for i in range(11):
    for ii in range(11):
        data[i][ii] = mapp(i,ii)



#print(data)


"""
par =  np.array(par)
std =  np.array(std)
corr = np.array(corr)

M, N = corr.shape # lines M, columns N ... vector = line

bb = [[1,2,3,4],
      [5,6,7,8],
      [9,0,1,2],
      [3,4,5,6],
      [7,8,9,0]]

aaa = np.array(bb).shape

print(aaa)

"""




"""
n = 200
#n=10
st = str("")
for i in range(1, n + 1):
    if i % 100 == 0: print(i)
    
    val2 = 0.2 + rn() * 0.4
    st +=  "aa{}".format(i)  + "   " + str(round(val2, 3)) + "   "
        
    for ii in range(1, i + 1):
        
        val = 0.1 + rn() * 0.5
        val2 = 0.2 + rn() * 0.4
        
        st +=  str(round(val, 3)) + "   "
    st += "1.0\n"




path = r"D:\_Trabalho\_Publicacoes_nossas\sampling\tests\corr.out"

f = open(path , "w")
f.write(st)
f.close()
"""



































