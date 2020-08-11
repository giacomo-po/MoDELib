# remember to "chmod +x plotStressStrain.py", then run as ./plotStressStrain.py
#!/usr/bin/python
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_table('F_0.txt', delim_whitespace=True,header=None)

plt.figure(0)
plt.plot(data[3],data[4])
plt.show()
