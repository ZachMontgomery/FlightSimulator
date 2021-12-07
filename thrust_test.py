import numpy as np

V = np.linspace(0.,300.,1000)

T0 = 1200.
T1 = -3.5
T2 = -.0022

T = T0 + T1 * V + T2 * V ** 2.

import matplotlib.pyplot as plt

plt.plot(V,T)
plt.show()
