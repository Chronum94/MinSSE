


class MinSSE:
  def __init__(lattice, beta, dilution, *, command = "./a.out", bins=100, msteps = 1000, isteps = 1000):
    self.lattice = lattice
    self.beta = beta
    self.dilution = dilution
    self.command = command
    self.bins = bins
    self.msteps = msteps
    self.isteps = isteps

  