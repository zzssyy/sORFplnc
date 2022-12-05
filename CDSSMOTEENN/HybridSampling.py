from CDSSMOTEENN import CDSSMOTE
import numpy as np
import random


def HybridSampling(X, y):
    print('Ours CDFSMOTE-WENN...')
    X_over, y_over, u = CDSSMOTE.CDSSMOTE(X, y, nomal=2, k=2, v=0.1, alpha=0.6, N=[0, 0], distance='cosin', rim_flag=2)
    enn = EditedNearestNeighbours(n_neighbors=5, kind_sel='mode')
    X_under, y_under = enn.fit_resample(X_over, y_over)

    return X_under, y_under, u