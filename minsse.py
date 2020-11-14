import os
import subprocess


class MinSSE:
    def __init__(self, 
        lattice,
        beta,
        dilution,
        *,
        directory=".",
        command="./a.out",
        bins=100,
        msteps=1000,
        isteps=1000,
    ):
        self.lattice = lattice
        self.beta = beta
        self.dilution = dilution

        self.directory = directory
        self.command = command

        self.bins = bins
        self.msteps = msteps
        self.isteps = isteps

    def run(self):
        if self.directory == '.':
          self.directory = "sse_rundir"
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        callstring = f"{self.command} {self.lattice[0]} {self.lattice[1]} {self.lattice[2]} {self.beta} {self.dilution} {self.bins} {self.msteps} {self.isteps}"
        subprocess.call(callstring, shell=True, cwd=self.directory)



callcommand = "/home/chronum/projects/SSE/MinSSE/a.out"
calc = MinSSE([6, 6, 6], 6.0, 0.0, command=callcommand, directory="Run1_realization1")

calc.run()