import numpy as np
import matplotlib.pyplot as plt

def Blahut_Arimoto(P, iteration_number=10000):
    R = []
    Q = []
    C = []
    C_output = []
    for i in range(0, iteration_number + 1):
        if i == 0:
            while True:
                random_r = np.random.rand(P.shape[1])
                if np.min(random_r) > 0:
                    R.append(random_r / sum(random_r))
                    break
        r = R[-1]
        q = np.zeros((P.shape[1], P.shape[0]))
        for y in range(P.shape[0]):
            s = np.dot(r, P[y, :].T)
            for x in range(P.shape[1]):
                q[x, y] = r[x] * P[y, x] / s
        Q.append(q)
        r = np.zeros(P.shape[1])
        for x in range(P.shape[1]):
            p = 1
            for y in range(P.shape[0]):
                if P[y, x] != 0:
                    p *= np.power(q[x, y], P[y, x])
            r[x] = p
        r /= sum(r)
        R.append(r)
        c = 0
        for x in range(P.shape[1]):
            for y in range(P.shape[0]):
                if P[y, x] != 0:
                    c += r[x] * P[y, x] * np.log2(q[x, y] / r[x])
        C.append(c)
        if i % 100 == 0:
            C_output.append(c)
    return r, c, C_output

iteration_number = 10000
P = np.matrix([[0.5, 0.7, 0.6], [0.3, 0.1, 0.05], [0.2, 0.2, 0.35]])
r, c, C_output = Blahut_Arimoto(P, iteration_number)
print(r, c)
fig = plt.figure(figsize=(10, 8))
plt.plot(np.array(list(range(0, iteration_number + 1, 100))), C_output)
plt.xlabel('iteration')
plt.ylabel('reachable capacity')
plt.show()