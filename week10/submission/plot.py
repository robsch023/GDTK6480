import re
import numpy as np
import matplotlib.pyplot as plt

fname = r"C:\Users\rober\OneDrive\Documents\GitHub\GDTK6480\week10\submission\wall-pressure.dat"
p_r = 500.0  # reservoir pressure, Pa

with open(fname) as f:
    lines = [ln.strip() for ln in f if ln.strip()]

# Header: first line, with '#', quotes and "N:" index prefixes stripped
header = lines[0].lstrip("#").replace('"', "").split()
header = [re.sub(r"^\d+:", "", name) for name in header]

data = np.array([[float(v) for v in ln.split()]
                 for ln in lines[1:] if not ln.startswith("#")])

x = data[:, header.index("pos.x")]
p = data[:, header.index("p")]

order = np.argsort(x)          # make sure x increases along the nozzle
x, p = x[order], p[order]

plt.figure(figsize=(8, 4.5))
plt.plot(x, p / p_r)
plt.xlabel("x (m)")
plt.ylabel("p / p$_r$")
plt.title("Wall pressure distribution, t = 5 ms")
plt.grid(True)
plt.tight_layout()
plt.savefig("wall-pressure.png", dpi=200)
plt.show()