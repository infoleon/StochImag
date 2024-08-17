# -*- coding: utf-8 -*-

"""

Read image as CSV and calculates the correlation (auto correlation) for the sampling:
    Correlation (based on auto-correlation)
    Parameters (already in possession)
    Std. Deviation (already in possession, but can be calculated based on the std. deviation)


"""


import pandas  as pd
import numpy as np

from copy import deepcopy

import sys
from operator import itemgetter


## INPUTS --------------------------------------------------------------

# kern = weight matrix, square matrix with ODD size (column and row)
kern = np.array([[1,1,1],
                 [1,1,1],
                 [1,1,1]])





# distance between each square pixel center
unit = 5.

# Inout files   -----------------------------
path = r"D:\work" + "\\"  # Folder to input files
file = r"csv_test_data.csv"
ful = path + file

# Output files -------------------------------
save_path  = r"D:\work" + "\\"
save_file0 = r"out_corr.csv"
save_file1 = r"out_par.csv"
save_file2 = r"out_std.csv"
save_file3 = r"out_param_coordinates.csv"
save_f0 = save_path + save_file0
save_f1 = save_path + save_file1
save_f2 = save_path + save_file2
save_f3 = save_path + save_file3

# Float value to be used where there is no data.
subst = -999

data = pd.read_csv(ful)
alles = data[["X", "Y", "Sand_00-30"]].values.tolist()
x = data["X"].values.tolist()
y = data["Y"].values.tolist()
s = data["Sand_60-90"].values.tolist()
ind_list = [i for i in range(len(s))]  # Labeling each cell, starting at 0

if len(x) != len(y) != len(s):
    print("fix input data, different quantity")
    sys.exit()

print("Inputs prepared\n")

#END OF INPUTS






## ---------------------------------------------------------------------
# FUNCTIONS
#%%&

def ind(kern = kern):
    """
    Returns a list with the "path" and the weight, for each value in the kern
    
    The path is the coordinate to "reach" to reach the value in the original matrix, which is needed for calculation.
    It considers zero coordinate as the current index.    
    """
    M = len(kern)
    start = int(-(M - 1)/2)  # path to the starting point
    res = [ [[0,0],0] for i in range(M * M)]  # Structure of the output
    
    indx = int(0) # general index
    #   i = y   ;   ii = x
    for i in range(M):
        for ii in range(M):
            yyy = start + i
            xxx = start + ii
            res[indx][0][0] = yyy
            res[indx][0][1] = xxx
            res[indx][1] = kern[i, ii]
            indx += int(1)
    return res

#END OF FUNCTIONS
#%%&
## ---------------------------------------------------------------------







# how many surrounding layers of pixels? (all directions, including diagonals)
sq_rad = (len(kern) - 1) / 2
sq_rad = int( sq_rad ) * int(2)  # times 2 (kern over kern)

max_x, min_x, max_y, min_y = max(x), min(x), max(y), min(y)
px = int(round((max_x - min_x) / unit)) + int(1)                 # number of pixels in x direction
py = int(round((max_y - min_y) / unit)) + int(1)
mat = np.zeros((py + int(2) * sq_rad, px + int(2) * sq_rad))     # Extended "matrix" filled with zeros, times 2 because it is 2 sides
mat[mat == 0.0] = subst
aver_m  = deepcopy(mat)
desv_m  = deepcopy(mat)
ind_mat = deepcopy(mat).astype(int)
raw_coord = deepcopy(mat).tolist()  # list with raw coordinates
# "mat" is the "picture" in the form of matrix

# Putting the data inside the created extended matrix
for i in range(len(x)):
    coor_x = int(round((x[i] - min_x) / unit))
    coor_y = int(round((y[i] - min_y) / unit))
    mat[coor_y + sq_rad][coor_x + sq_rad] = s[i]
    ind_mat[coor_y + sq_rad][coor_x + sq_rad] = int(ind_list[i]) # Label names
    raw_coord[coor_y + sq_rad][coor_x + sq_rad] = [x[i], y[i]]  # X, Y

# path and val, output
# [[-1, -1], 1]





# Creating correlation matrix
corr = []
#kern_t = deepcopy(kern)  # temporary kern for 
for i in range(len(mat)):
    for ii in range(len(mat[i])):
        #   i = y   ;   ii = x
        #  i, ii  ->  Coord base cell
        if mat[i][ii] != subst:
            for iii in ind(kern):  # -------------------------------------------------------------------------------------
                #  y, x  ->  Coord 1 kern cell
                y = i  + iii[0][0]
                x = ii + iii[0][1]
                if mat[y][x] != subst:
                    
                    # -------------------------------------------------------------------------------------
                    # AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE
                    aver1 = 0.0
                    aver2 = 0.0
                    n_ave = int(0)
                    for iv in ind(kern):  # Enter kern ------------------------------------------------
                        # y_b, x_b -> both coordinates
                        y_b1 = i  + iv[0][0]
                        x_b1 = ii + iv[0][1]
                        
                        y_b2 = y  + iv[0][0]
                        x_b2 = x  + iv[0][1]
                        
                        if (mat[y_b1][x_b1] != subst) and (mat[y_b2][x_b2] != subst):
                            
                            val1 = mat[y_b1][x_b1] # value of the property
                            val2 = mat[y_b2][x_b2]
                            
                            aver1 += val1
                            aver2 += val2
                            
                            n_ave += int(1)
                    
                    if n_ave > int(0):
                        aver1 = aver1 / n_ave # averages
                        aver2 = aver2 / n_ave
                            
                    # -------------------------------------------------------------------------------------
                    # DEVIATION DEVIATION DEVIATION DEVIATION DEVIATION DEVIATION DEVIATION DEVIATION
                    desv1 = 0.0
                    desv2 = 0.0
                    n_desv = int(0)
                    for iv in ind(kern):  # Enter kern ------------------------------------------------
                        # y_b, x_b -> both coordinates
                        y_b1 = i  + iv[0][0]
                        x_b1 = ii + iv[0][1]
                        
                        y_b2 = y  + iv[0][0]
                        x_b2 = x  + iv[0][1]
                        
                        if (mat[y_b1][x_b1] != subst) and (mat[y_b2][x_b2] != subst):
                            
                            val1 = mat[y_b1][x_b1]  # value of the property
                            val2 = mat[y_b2][x_b2]
                            
                            desv1 += (val1 - aver1) ** 2
                            desv2 += (val2 - aver2) ** 2
                            
                            n_desv += int(1)
                      
                    if n_desv > int(0):
                        desv1 = (desv1 / n_desv) ** 0.5 # deviaitons
                        desv2 = (desv2 / n_desv) ** 0.5
                    
                    # -------------------------------------------------------------------------------------
                    # CORRELATION CORRELATION CORRELATION CORRELATION CORRELATION CORRELATION CORRELATION
                    _sum  = 0.0  # sum for correlation
                    sum_n = int(0)  # sum of values used
                    
                    ind1  = ind_mat[i][ii] # cell name (label)
                    ind2  = ind_mat[y][x]  # cell name (label)
                    
                    for iv in ind(kern): # Enter kern ------------------------------------------------
                        # y_b, x_b -> both coordinates
                        y_b1 = i  + iv[0][0]
                        x_b1 = ii + iv[0][1]
                        
                        y_b2 = y  + iv[0][0]
                        x_b2 = x  + iv[0][1]
                        if (mat[y_b1][x_b1] != subst) and (mat[y_b2][x_b2] != subst):
                            
                            val1 = mat[y_b1][x_b1] # value of the property inside kernel
                            val2 = mat[y_b2][x_b2]
                            
                            _sum += (val1 - aver1) * (val2 - aver2)
                            sum_n += int(1)
                    
                    # BU data
                    val    = mat[i][ii] # ind1 value
                    _coord = raw_coord[i][ii]
                    
                    if sum_n > int(0):
                        
                        _corr = _sum / (sum_n * desv1 * desv2)
                        
                        # The coordinate is relative to the first value (ind1)
                        #          0    1           2            3      4      5     6     7      
                        _temp = [ind1, ind2, round(_corr, 4), _coord, desv1 , val, aver1, sum_n]
                        corr.append(_temp)









_ind_t = []
for i in range(len(corr)):
    _ind_t.append(corr[i][0])

_ind = sorted(set(_ind_t))  # ind sorted
N = len(_ind)

_l = deepcopy(corr)
sorted(_l, key=itemgetter(0))  # sorted corr

_data_out = []     # list to build the parameter and std deviation output
for i in _ind:
    for ii in _l:
        if i == ii[0]:
            _data_out.append(ii)
            break

corr_mat = np.zeros((N,N))
for i in range(len(corr)):
    _ind1 = corr[i][0]
    _ind2 = corr[i][1]
    corr_mat[_ind1][_ind2] = corr[i][2]

text = ""
for i in range(len(corr_mat)):
    if (i == len(corr_mat) - int(1)):
        text += "par" + str(i)
    else:
        text += "par" + str(i) + ", "
par_str = deepcopy(text) + "\n"
std_str = deepcopy(text) + "\n"

for i in range(len(corr_mat)):
    text += "\n"
    for ii in range(len(corr_mat[i])):
        if ii == (len(corr_mat[i])-1):
            text += str(corr_mat[i][ii])
        else:
            text += str(corr_mat[i][ii]) + ", "

bu_ccord = []
for i in range(N):
    par_str += str(_data_out[i][5]) + ", "
    std_str += str(_data_out[i][4]) + ", "
    bu_ccord.append([ _data_out[i][0] , _data_out[i][3][0], _data_out[i][3][1] ])

par_str = par_str[:-2]
std_str = std_str[:-2]


with open(save_f0, "w+") as f:
    f.write(text)

with open(save_f1, "w+") as f:
    f.write(par_str)

with open(save_f2, "w+") as f:
    f.write(std_str)

    

df_bu = pd.DataFrame(bu_ccord, columns=[["param","X","Y"]])
df_bu.to_csv(save_f3, index = False)























