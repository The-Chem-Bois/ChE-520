import numpy as np
from Hmw1 import *

coefficients = np.array([[16.1156,3395.57,-59.44],[16.1156,3395.57,-59.44],[16.1156,3395.57,-59.44]]) 
print(antoine(361.10,coefficients, output="P"))