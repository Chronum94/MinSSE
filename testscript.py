import subprocess
import os



if not os.path.exists("TrialDir"):
  os.makedirs("TrialDir")

subprocess.call("/home/chronum/projects/SSE/MinSSE/a.out 10 10 10 16.0 0.3 10 1000 1000",
shell=True, cwd="TrialDir/")