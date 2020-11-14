import numpy as np
import subprocess
import os


def run_single_beta(i, beta):
  dirname = "Test2D_16x16_highbetasweep_{}_dil0.40".format(i)
  if not os.path.exists(dirname):
    os.makedirs(dirname)

  subprocess.call("/home/chronum/projects/SSE/MinSSE/a.out 16 16 1 {} 0.40 10000 5000 5000".format(beta),
  shell=True, cwd=dirname)


beta = np.logspace(1, 4, num=6, base=2.0)

from multiprocessing import Pool

with Pool(6) as p:
  p.starmap(run_single_beta, [(i, b) for i, b in enumerate(beta)])
