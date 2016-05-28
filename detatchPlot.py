# -*- coding: utf-8 -*-
"""
Created on Tue May 24 18:42:03 2016

@author: Mark
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]
data = np.load(filename)
plt.plot(data['x'], data['y'])
plt.show()