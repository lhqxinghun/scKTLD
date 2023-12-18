# encoding=utf8
import sys
from sklearn import manifold
from sklearn import datasets
import numpy as np

import scKTLD.cpd_utils as cpd
#import ticc


# ticc <Multivariate Time Series data  clustering>

def cpDector (signal,  method= "KernelCPD", min_size = 3, **kwargs):
    print("Starting TAD boundary detection")
    algo = None
    if method == "KernelCPD":
        kernel = "rbf"
        min_size = 3
        if "kernel" in kwargs:
            kernel = kwargs["kernel"]
        if "min_size" in kwargs:
            min_size = kwargs["min_size"]
        algo = cpd.KernelCPD(kernel=kernel, jump=1, min_size=min_size).fit(signal)

    if method == "Pelt":
        model = "rbf"
        min_size = 3
        if "model" in kwargs:
            model = kwargs["model"]
        if "min_size" in kwargs:
            min_size = kwargs["min_size"]
        algo = cpd.Pelt(model=model, jump=1, min_size=min_size).fit(signal)

    if method == "Binseg":
        model = "rbf"
        min_size = 3
        if "model" in kwargs:
            model = kwargs["model"]
        if "min_size" in kwargs:
            min_size = kwargs["min_size"]
        algo = cpd.Binseg(model=model, jump=1, min_size=min_size).fit(signal)

    if method == "BottomUp":
        model = "rbf"
        min_size = 3
        if "model" in kwargs:
            model = kwargs["model"]
        if "min_size" in kwargs:
            min_size = kwargs["min_size"]
        algo = cpd.BottomUp(model=model, jump=1, min_size=min_size).fit(signal)

    if algo == None:
        print("Error: Unkown search method")
        sys.exit(-1)
    else:
        pen = 2
        if "pen" in kwargs:
            pen = kwargs["pen"]
        result = algo.predict(pen = pen)
        result = [0]+result
        result[-1] = result[-1]-1
    print("TAD boundary detection Done!")

    return np.array(result)


        

