# scKTLD

## 1. Introduction
scKTLD is a method designed for the identification of TAD-like domains on single-cell Hi-C data. It treats the Hi-C contact matrix as a graph, embeds its structures into a low-dimensional space by combining sparse matrix factorization and spectral propagation, and identifies the TAD-like domains in the embedding space via a kernel model optimized by PELT
![image](https://github.com/lhqxinghun/scKTLD/blob/main/data/exp-sc/overview.PNG)

## 2. Installation & Example

**2.1 OS**
- ubuntu 18.04

**2.2 Required Python Packages**
Make sure that all the packages listed in the *requirements.txt* are installed.

- Python >= 3.6
- scipy >= 1.5.2
- numpy >= 1.18.0
- ticc >= 0.1.6
- networkx >= 2.5
- scikit-learn >= 0.24.2

**2.3 Install from Github**

(1) Download the folder *scKTLD* by git clone
```
$ git clone https://github.com/lhqxinghun/scKTLD/
```
(2) Install the package *scKTLD* with the following command:
```
$ conda create -n scKTLD python=3.6
$ pip install Cython
$ cd scKTLD
$ pip install . #or you can try python setup.py install 
```

**2.4 Run example**
```
$ cd scKTLD
$ python example.py
# If it works properly, the following figure will be shown
```
![image](https://github.com/lhqxinghun/scKTLD/blob/main/data/exp-sc/Result.png)

More detailed examples can be find in the jupyter notebook *example.ipynb*

## 3. Usage
### 3.1 Key function *callTLD*
#### Input parameters:
- **graph** np.ndarray, the dense format of a contact matrix.
- **dimension** int, dimension of the embedding vectors of nodes.
- **penalty** float, penalty constant during changepoint detection.
- **brecon** bool, whether or not to return the reconstructed Hi-C map.
## 4. Contact
hongqianglv@mail.xjtu.edu.cn OR liuerhu@stu.xjtu.edu.cn
