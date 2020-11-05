import numpy as np

beta = 2
nbonds = 2 * 16

L = 10

opstring = np.zeros(L)

spinpar = np.random.choice([0, 1], size=nbonds)

n = 0

for p in range(L):
    paccept = 0.5 * beta * nbonds
    pdelete = 1 / paccept
    if opstring[p] == 0:
        b = np.random.randint(0, nbonds)
        if not spinpar[b]:
            if np.random.rand() * (L - n) <= paccept:
                opstring[p] = 2 * b
                n += 1
    elif opstring[p] % 2 == 0:
        if np.random.rand() * (L - n + 1) <= pdelete:
            opstring[p] = 0
            n -= 1
    

print(len(opstring), len(np.nonzero(opstring)[0]))