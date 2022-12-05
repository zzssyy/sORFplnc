import numpy as np
from sklearn import metrics
import scipy.stats as stats
from skfeature.function.similarity_based import fisher_score
import networkx as nx
import copy
import random as ra
from CDSSMOTEENN import RIMquantifier
from scipy import spatial
from CDSSMOTEENN.IntraClassImbalance import IntraClassImbalance
from sklearn.feature_selection import VarianceThreshold
from LaplacianScoreforFeatureSelection.lp_score import *

def get_AdjacencyMatrix(matrix):
    # Construct an adjacency matrix A
    results = np.zeros(shape=(len(matrix), len(matrix)), dtype='float64')
    for i in range(0, len(matrix)-1):
        if i == 0:
            continue
            # results[matrix[i]][matrix[i+1]] = 1
        else:
            results[matrix[i]][matrix[i-1]] = 1
            # results[matrix[i]][matrix[i+1]] = 1
    return results

def FeatureRanking(X, y):
    X_T = X.T

    # laplacian score
    LS = LaplacianScore(X, neighbour_size=16, t_param=2)
    LS = feature_ranking(LS)

    # VarianceThreshold
    vt = VarianceThreshold()
    vt.fit(X)
    feaures = [i for i in range(len(X.T))]
    VTs = [(x, y) for x, y in zip(feaures, vt.variances_)]
    VTs = np.array([i[0] for i in sorted(VTs, key=lambda x: abs(x[1]), reverse=True)])

    # Fisher Score
    FS = fisher_score.fisher_score(X, y, mode='rank')

    # Mutual Information
    MI = list()
    for i in X_T:
        MI_result = metrics.normalized_mutual_info_score(i, y)
        MI.append(MI_result)
    MI = np.argsort(np.array(MI))

    # Correlation Score
    r = list()
    for i in X_T:
        results = list()
        for j in X_T:
            result = stats.pearsonr(i, j)
            results.append(result[0])
        r.append(results)
    CS = np.sum(r, axis=0)
    CS = np.argsort(CS)

    # Eigenvector Centrality
    A1 = get_AdjacencyMatrix(LS)
    A2 = get_AdjacencyMatrix(VTs)
    A3 = get_AdjacencyMatrix(FS)
    A4 = get_AdjacencyMatrix(MI)
    A5 = get_AdjacencyMatrix(CS)
    value = np.array(A1) + np.array(A2) + np.array(A3) + np.array(A4) + np.array(A5)
    vertices_s = [i for i in range(0, len(value))]
    vertices_e = [i for i in range(0, len(value))]
    A = nx.DiGraph()
    for i in range(np.size(vertices_s)):
        for j in range(np.size(vertices_e)):
            A.add_weighted_edges_from([(vertices_s[i], vertices_e[j], value[i, j])])
    eigenvector = nx.pagerank(A, max_iter=1000)
    eigenvector = sorted(eigenvector.items(), key=lambda d: d[1], reverse=True)
    order = []
    for key in eigenvector:
        order.append(key[0])

    return order

def get_Distance(a ,b, nomal, distance='Minkowski'):
    if distance == 'Minkowski':
        #Minkowski distance
        MD = np.linalg.norm((a-b), nomal)
        return MD
    else:
        # cosin distance
        CS = 1 - spatial.distance.cosine(a, b)
        return CS

def get_TopK(matrix, k):
    t = copy.deepcopy(matrix)
    max_number = []
    max_index = []
    for _ in range(k):
        number = min(t)
        index = t.index(number)
        a = float('inf')
        t[index] = a
        max_number.append(number)
        max_index.append(index)
    t = []
    return max_index

def IMOWADSMOTE(X, weight, nomal, k, N, v, distance='cosin'):
    M = list()
    MM = list()
    for i in X:
        result = i * weight
        M.append(result)
        MM.append(i.tolist())
    for i in range(0, len(M)):
        IMOWAD = list()
        for j in range(0, len(M)):
            if i != j:
                distance = get_Distance(M[i], M[j], nomal, distance='cosin')
                IMOWAD.append(distance)
            else:
                distance = float('inf')
                IMOWAD.append(distance)
        index = get_TopK(IMOWAD, k)
        T = [M[i] for i in index]
        MM = synthesissample(T, v, X[i], MM, N[i])
    return MM

def synthesissample(T, v, xi, MM, N):
    ra.shuffle(T)
    for k in range(0, N):
        xk = T[0]
        xi = np.array(xi)
        xk1 = xi + v*(xk - xi)
        MM.append(xk1.tolist())
        T = np.delete(np.array(T), 0)
    return MM

def feature_selection(M, fs, u):
    z = M.T
    x = np.array([z[i] for i in range(len(z)) if i in u[:fs]]).T
    return x

def FSSMOTE(X, y, nomal, k, v, alpha, N, distance, rim_flag):
    X, y = np.array(X), np.array(y)
    X = X[np.where(y == 0)].tolist() + X[np.where(y == 1)].tolist()
    y = y[np.where(y == 0)].tolist() + y[np.where(y == 1)].tolist()
    X, y = np.array(X), np.array(y)
    # Feature ranking methods
    u = FeatureRanking(X, y)
    fs = int(0.5 * len(u))
    X = feature_selection(X, fs, u)
    weight = RIMquantifier.RIM(u[:fs], alpha, n=len(X.T), rim_flag=2)
    dense = IntraClassImbalance(X, y, n_clusters=4)
    N = np.array(dense)

    # Judge minority class and amount of oversampling
    if len(y[y == 1]) != len(y[y == 0]):
        if len(y[y == 1]) > len(y[y == 0]):
            X_min, X_maj = X[np.where(y == 0)], X[np.where(y == 1)]
            y_min, y_maj = y[np.where(y == 0)], y[np.where(y == 1)]
            N_min = N[np.where(y == 0)]
            # N = int((len(y[y == 1]) - len(y[y == 0])) / len(y[y == 0]))
        elif len(y[y == 1]) < len(y[y == 0]):
            X_min, X_maj = X[np.where(y == 1)], X[np.where(y == 0)]
            y_min, y_maj = y[np.where(y == 1)], y[np.where(y == 0)]
            N_min = N[np.where(y == 1)]
            # N = int((len(y[y == 0]) - len(y[y == 1])) / len(y[y == 1]))

        # Core code of FS-SMOTE
        M = IMOWADSMOTE(X_min, weight, nomal, k, N_min, v, distance='cosin')
        M_len = np.array(list(set(y_min)) * len(M))
        X = np.row_stack([M, X_maj])
        y = np.concatenate([M_len, y_maj])
        return X, y, u[:fs]
    else:
        print('This is not an imbalanced dataset problem. Only selecting feature.')
        return X, y, u[:fs]

# if __name__ == "__main__":
#     FSSMOTE(X=[[1, 2, 3, 4, 5, 1], [0, 0, 1, 1, 2, 1], [2, 2, 3, 4, 5, 7], [3, 0, 1, 1, 2, 4], [1, 2, 3, 4, 2, 7],
#                    [1, 0, 1, 1, 2, 7], [4, 5, 6, 9, 7, 10]], y=[0, 0, 0, 0, 0, 1, 1], nomal=2, k=2, v=0.1, alpha=0.6)