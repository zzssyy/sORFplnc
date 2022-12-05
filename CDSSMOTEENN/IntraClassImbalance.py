from sklearn import preprocessing
from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
import numpy as np
from scipy import spatial
from CDSSMOTEENN import CDSSMOTE

def distance(a, b):
    dis = 1 - spatial.distance.cosine(a, b)
    return dis

def average_distance(lst):
    result = []
    sum = 0
    count = 0
    for i in lst:
        sum = sum + i[1]
        count = count + 1
    sum = sum / count
    for z in lst:
        index = z[0]
        fre = sum / z[1]
        s = (index, fre)
        result.append(s)
    return result

def IntraClassImbalance(X, y, n_clusters=2):
    kmeans = KMeans(init='k-means++', n_clusters=n_clusters, random_state=10).fit(preprocessing.normalize(X))
    k_labels = kmeans.predict(X)
    flag = [i for i in range(len(y)-1) if y[i] != y[i+1]][0]
    center_ps = kmeans.cluster_centers_
    set_k_labels = set(k_labels)
    index = []
    for i in set_k_labels:
        num = list(np.where(k_labels == i)[0])
        index.append(num)
    count = []
    for i in index:
        posC = 0
        negC = 0
        for j in i:
            if j <= flag:
                negC = negC + 1
            else:
                posC = posC + 1
        count.append([negC, posC])

    count_min = []
    for i in count:
        if 0.15 * i[0] >= i[1]:
            count_min.append(0)
        else:
            count_min.append(i[1]/len(X[y == 1]))

    dis = []
    for j in range(0, len(center_ps)):
        result = []
        for i in range(0, len(k_labels)):
            if k_labels[i] == j and i > flag:
                result.append((i, distance(X[i], center_ps[j])))
            else:
                continue
        dis.append(result)

    fre = []
    for i in dis:
        result = []
        if len(i) == 0:
            result = result + []
        else:
            result = result + average_distance(i)
        fre.append(result)

    dense = [0]*len(X)
    gap = 2*(flag+1)-len(X)
    for i in range(0, len(count_min)):
        if len(fre[i]) == 0:
            continue
        else:
            for j in fre[i]:
                if gap*j[1]//(len(X)-flag-1) == 0:
                    m = 1
                else:
                    m = gap*j[1]//(len(X)-flag-1)
                dense[j[0]] = int(m)
    return dense

# if __name__ == "__main__":
    # iris = load_iris()
    # x = iris['data'][:70]
    # y = iris['target'][:70]
    # dense = IntraClassImbalance(x, y, n_clusters=4)
    # x, y = CDSSMOTEENN.CDSSMOTE(x, y, nomal=2, k=2, v=0.1, alpha=0.6, N=np.array(dense))
    # print(x)
    # print(y)
    # print(len(x))
    # print(len(y))