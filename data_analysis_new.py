import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.linear_model import LinearRegression


run_regression = True

phi_array = np.genfromtxt("simulation_phi.csv", delimiter =",")
rphi_array = np.genfromtxt("simulation_rphi.csv", delimiter =",")
dz_rphi_array = np.genfromtxt("simulation_dz_rphi.csv", delimiter =",")
ddz_rphi_array = np.genfromtxt("simulation_ddz_rphi.csv", delimiter =",")
dz_phi_array = np.genfromtxt("simulation_dz_phi.csv", delimiter =",")
ds_rphi_array = np.genfromtxt("simulation_ds_rphi.csv", delimiter =",")
r_dz_phi_array = dz_phi_array - .5*phi_array



max_dz_rphi_array = np.nanmax(np.abs(dz_rphi_array),axis=1)
max_ds_rphi_array = np.nanmax(np.abs(ds_rphi_array),axis=1)


ns = phi_array.shape[0]
nz = phi_array.shape[1]
print(ns)
rangez = 1


ds = .0001 ## Needs to be changed as simulation parameters are updated
ranges = 10  ## Needs to be changed as simulation parameters are updated

z_array = np.linspace(-1,0, nz)

s_array = np.linspace(0, ranges, ns)


if(run_regression):
    reg = LinearRegression().fit(s_array[400:ns].reshape(-1, 1), np.log(max_dz_rphi_array[400:ns]))
    print(reg.coef_)
    reg1 = LinearRegression().fit(s_array[400:ns].reshape(-1, 1), np.log(max_ds_rphi_array[400:ns]))
    print(reg1.coef_)


plt.plot(s_array,np.log(max_dz_rphi_array))
plt.waitforbuttonpress()
plt.clf()

for i in range(ns):

    plt.plot(z_array, rphi_array[i,:])
    plt.plot(z_array, dz_rphi_array[i,:])
    plt.plot(z_array,ds_rphi_array[i,:])
    plt.plot(z_array,dz_phi_array[i,:])

    plt.legend(['rphi','dz_rphi','ds_rphi','dz_phi'])
    plt.title(s_array[i])
    plt.draw()
    plt.waitforbuttonpress()
    plt.clf()
