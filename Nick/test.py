"""
This is the testing script
"""
import numpy as np
import matplotlib.pyplot as plt
from fields import SHField


import fields as f

g = np.array([[10000, 10000]])
h = np.array([[0, -10000]])
field1 = SHField(1000000, g, h, 0, 0)

field1.plotField(-100000, 100000, -100000, 100000, -100000, 100000)


plt.show()