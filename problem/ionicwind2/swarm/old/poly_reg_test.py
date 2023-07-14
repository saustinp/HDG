import numpy as np
import matplotlib.pyplot as plt
from scipy import odr

# Simple script to test the capability of the python regression functions
x = np.linspace(0.0, 5.0)
y = np.sin(x)
poly_model = odr.polynomial(1)  # using third order polynomial model
data = odr.Data(x, y)
odr_obj = odr.ODR(data, poly_model)
output = odr_obj.run()  # running ODR fitting
poly = np.poly1d(output.beta[::-1])
poly_y = poly(x)
plt.plot(x, y, label="input data")
plt.plot(x, poly_y, label="polynomial ODR")
plt.legend()
plt.show()