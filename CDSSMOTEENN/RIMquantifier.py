import numpy as np

def RIM(u, alpha, n, rim_flag):

    u1 = [(i + 1) / n for i in u]
    u2 = [i / n for i in u]

    if rim_flag == 0:
        # Basic RIM quantifier
        b1 = np.power(u1, alpha)
        b2 = np.power(u2, alpha)
        BRIM = [i-j for i, j in zip(b1, b2)]
        # print(BRIM)
        return BRIM

    if rim_flag == 1:
        # Quadratic RIM quantifier
        q1 = [1/(1-alpha*i) for i in np.power(u1, 0.5)]
        q2 = [1/(1-alpha*i) for i in np.power(u2, 0.5)]
        QRIM = [i-j for i, j in zip(q1, q2)]
        # print(QRIM)
        return QRIM

    if rim_flag == 2:
        # Exponential RIM quantifier
        e1 = np.exp(-alpha * np.array(u1, dtype=float))
        e2 = np.exp(-alpha * np.array(u2, dtype=float))
        ERIM = [i - j for i, j in zip(e1, e2)]
        # print(ERIM)
        return ERIM

    if rim_flag == 3:
        # Trigonometric RIM quantifier
        t1 = np.arcsin(alpha * np.array(u1, dtype=float))
        t2 = np.arcsin(alpha * np.array(u2, dtype=float))
        TRIM = [i - j for i, j in zip(t1, t2)]
        # print(TRIM)
        return TRIM

    if rim_flag == 4:
        ones = [1] * len(u)
        return ones