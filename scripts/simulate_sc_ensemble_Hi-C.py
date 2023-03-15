import numpy as np
import pandas as pd
import os

path_input = "/media/biology/datadisk/liuerhu/scTAD/simulatedata/K562-070/outputs_test_800/chr8_10000_10499/models"
xyzname = [f for f in os.listdir(path_input) if f.endswith('.xyz')]
model_number = [int(str.split(f, '.')[1]) for f in xyzname]

for mn in model_number:
    ##load 10kb resolution 3D model from IMP
    xyz=np.array(pd.read_csv(path_input+"/model."+ str(mn) + ".xyz",delimiter='\t',skiprows=0,header=None).loc[:,[2,3,4]])

    ##generate distance matrix
    dist_simu=np.zeros([np.shape(xyz)[0],np.shape(xyz)[0]])
    for i in range(np.shape(xyz)[0]):
        for j in range(np.shape(xyz)[0]):
            dist_simu[i,j]=((xyz[i,0]-xyz[j,0])**2+(xyz[i,1]-xyz[j,1])**2+(xyz[i,2]-xyz[j,2])**2)**(0.5)

    ##Expect contacts number of the simulated scHi-C data
    cov=1000

    ##Only pairs of genome loci with a Euclidean distance less than a threshold (dis) in the model can have contacts
    dis=500
    ##simulate 40kb resolution scHi-C data
    size=np.shape(dist_simu)[0]//5
    C=np.zeros([size,size])
    A=dist_simu
    np.random.seed(0)
    weight=np.sum(np.triu((dis-A)*(A<dis),0))
    rate=cov/weight
    for i in range(size):
        for j in range(i+1,size):
            for x in range(5*i,5*i+5):
                for y in range(5*j,5*j+5):
                    if A[x,y]<dis and y>=x :
                        C[i,j]+=np.random.binomial(1,rate*(dis-A[x,y]))
    C=C+np.triu(C,0).T
    np.savetxt(path_input+"/thre_500/model."+str(mn)+".mat", C, fmt='%d')


    ##Only pairs of genome loci with a Euclidean distance less than a threshold (dis) in the model can have contacts
    dis=750
    ##simulate 40kb resolution scHi-C data
    size=np.shape(dist_simu)[0]//5
    C=np.zeros([size,size])
    A=dist_simu
    np.random.seed(0)
    weight=np.sum(np.triu((dis-A)*(A<dis),0))
    rate=cov/weight
    for i in range(size):
        for j in range(i+1,size):
            for x in range(5*i,5*i+5):
                for y in range(5*j,5*j+5):
                    if A[x,y]<dis and y>=x :
                        C[i,j]+=np.random.binomial(1,rate*(dis-A[x,y]))
    C=C+np.triu(C,0).T
    np.savetxt(path_input+"/thre_750/model."+str(mn)+".mat", C, fmt='%d')

    ##Only pairs of genome loci with a Euclidean distance less than a threshold (dis) in the model can have contacts
    dis=1000
    ##simulate 40kb resolution scHi-C data
    size=np.shape(dist_simu)[0]//5
    C=np.zeros([size,size])
    A=dist_simu
    np.random.seed(0)
    weight=np.sum(np.triu((dis-A)*(A<dis),0))
    rate=cov/weight
    for i in range(size):
        for j in range(i+1,size):
            for x in range(5*i,5*i+5):
                for y in range(5*j,5*j+5):
                    if A[x,y]<dis and y>=x :
                        C[i,j]+=np.random.binomial(1,rate*(dis-A[x,y]))
    C=C+np.triu(C,0).T
    np.savetxt(path_input+"/thre_1000/model."+str(mn)+".mat", C, fmt='%d')

    ##Expect contacts number of the simulated ensemble Hi-C data
    cov=350000

    ##simulate 50kb resolution ensemble Hi-C data
    C=np.zeros([size,size])
    A=dist_simu
    np.random.seed(0)
    rate=cov/np.sum(np.triu(1/dist_simu,1))
    for i in range(size):
        for j in range(i+1,size):
            for x in range(5*i,5*i+5):
                for y in range(5*j,5*j+5):
                    C[i,j]+=np.random.poisson(rate/A[x,y])
    C=C+np.triu(C,1).T
    np.savetxt(path_input+"/ensemble/model."+str(mn)+".mat", C, fmt='%d')

