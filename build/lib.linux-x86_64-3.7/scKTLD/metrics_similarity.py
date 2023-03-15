import numpy as np
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import adjusted_rand_score

def ami(C1, C2):
    return adjusted_mutual_info_score(C1, C2)

def ari(C1, C2):
    return adjusted_rand_score(C1, C2)

def moc(P, Q):
    NP = P.shape[0]
    NQ = Q.shape[0]

    if(NP==NQ & NP == 1):
        moc = 1
    else:
        moc = 0
        for i in range(NP):
            for j in range(NQ):
                Pi = np.arange(P[i][0], P[i][1])
                Qj = np.arange(Q[j][0], Q[j][1])
                Fij = np.intersect1d(Pi, Qj)
                moc = moc + (len(Fij*Fij)**2/(len(Pi)*len(Qj)))
        moc = (moc-1)/((NP*NQ)**(1/2)-1)
    
    return moc

def vi(X, Y):
    k = X.shape[0]
    l = Y.shape[0]

    n1 = 0
    for i in range(k):
        Xi = np.arange(X[i, 0],  X[i, 1])
        n1 = n1 + len(Xi)
    n2 = 0
    for j in range(l):
        Yj = np.arange(Y[j, 0],  Y[j, 1])
        n2 = n2 + len(Yj)

    if (n1 != n2):
        print("Warning: the number of elements in X and Y are inconsistent, using the former")
    n = n1

    vi = 0
    for i in range(k):
        for j in range(l):
            Xi = np.arange(X[i, 0],  X[i, 1])
            Yj = np.arange(Y[j, 0],  Y[j, 1])
            pi = len(Xi)/n
            qj = len(Yj)/n
            rij = len(np.intersect1d(Xi, Yj))/n
            if(rij!=0):
                vi = vi -rij*(np.log2(rij/pi) + np.log2(rij/qj))
    return vi


from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import adjusted_rand_score