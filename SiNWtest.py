import numpy as np
from SiNWfunction import bs

num_points = 10
kk = np.linspace(0, 0.57, num_points, endpoint=True)
Eg, bstruct = bs(path='c:\users\sammy\desktop\NanoNet\input_samples', kk=kk, flag=True)

print Eg
