
# encoding=utf8
import numpy as np
import networkx as nx
import scipy.sparse
import scipy.sparse as sp
from scipy import linalg
from scipy.special import iv

from sklearn import preprocessing
from sklearn.utils.extmath import randomized_svd

import time

import sklearn.metrics.pairwise as skmetrics

class ProNE():
	def __init__(self, graph, dimension):
		self.dimension = dimension

		
		'''self.G = nx.read_edgelist(self.graph, nodetype=int, create_using=nx.DiGraph())
		self.G = self.G.to_undirected()
		self.node_number = self.G.number_of_nodes()
		matrix0 = scipy.sparse.lil_matrix((self.node_number, self.node_number))

		for e in self.G.edges():
			if e[0] != e[1]:
				matrix0[e[0], e[1]] = 1
				matrix0[e[1], e[0]] = 1
		self.matrix0 = scipy.sparse.csr_matrix(matrix0)'''
		
		self.G=graph
		self.node_number=len(self.G)
		self.matrix0 = scipy.sparse.csr_matrix(self.G)
		print(self.G.shape)

	def get_embedding_rand(self, matrix):
		# Sparse randomized tSVD for fast embedding
		t1 = time.time()
		l = matrix.shape[0]
		smat = scipy.sparse.csc_matrix(matrix)  # convert to sparse CSC format
		#print('svd sparse', smat.data.shape[0] * 1.0 / l ** 2)
		U, Sigma, VT = randomized_svd(smat, n_components=self.dimension, n_iter=5, random_state=None)
		U = U * np.sqrt(Sigma)
		U = preprocessing.normalize(U, "l2")
		print('sparsesvd time', time.time() - t1)
		return U

	def get_embedding_dense(self, matrix, dimension):
		# get dense embedding via SVD
		t1 = time.time()
		U, s, Vh = linalg.svd(matrix, full_matrices=False, check_finite=False, overwrite_a=True)
		U = np.array(U)
		U = U[:, :dimension]
		s = s[:dimension]
		s = np.sqrt(s)
		U = U * s
		U = preprocessing.normalize(U, "l2")
		#print('densesvd time', time.time() - t1)
		print('Node embedding done!')
		return U

	def pre_factorization(self, tran, mask):
		# Network Embedding as Sparse Matrix Factorization
		t1 = time.time()
		l1 = 0.75
		C1 = preprocessing.normalize(tran, "l1")
		neg = np.array(C1.sum(axis=0))[0] ** l1

		neg = neg / neg.sum()

		neg = scipy.sparse.diags(neg, format="csr")
		neg = mask.dot(neg)
		print("neg", time.time() - t1)

		C1.data[C1.data <= 0] = 1
		neg.data[neg.data <= 0] = 1

		C1.data = np.log(C1.data)
		neg.data = np.log(neg.data)

		C1 -= neg
		F = C1
		features_matrix = self.get_embedding_rand(F)
		return features_matrix

	def chebyshev_gaussian(self, A, a, order=10, mu=0.5, s=0.5):
		# NE Enhancement via Spectral Propagation
		print('Chebyshev expansion')
		t1 = time.time()

		if order == 1:
			return a

		A = sp.eye(self.node_number) + A
		DA = preprocessing.normalize(A, norm='l1')
		L = sp.eye(self.node_number) - DA

		M = L - mu * sp.eye(self.node_number)

		Lx0 = a
		Lx1 = M.dot(a)
		Lx1 = 0.5 * M.dot(Lx1) - a

		conv = iv(0, s) * Lx0
		conv -= 2 * iv(1, s) * Lx1
		for i in range(2, order):
			Lx2 = M.dot(Lx1)
			Lx2 = (M.dot(Lx2) - 2 * Lx1) - Lx0
			#         Lx2 = 2*L.dot(Lx1) - Lx0
			if i % 2 == 0:
				conv += 2 * iv(i, s) * Lx2
			else:
				conv -= 2 * iv(i, s) * Lx2
			Lx0 = Lx1
			Lx1 = Lx2
			del Lx2
			#print('Bessell time', i, time.time() - t1)
		mm = A.dot(a - conv)
		emb = self.get_embedding_dense(mm, self.dimension)
		return emb