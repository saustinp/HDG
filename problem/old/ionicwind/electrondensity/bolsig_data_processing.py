import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import odr

def create_default_fig(title, xlabel, ylabel, figsize=(16,10), grid=1):
    fig, ax = plt.subplots(figsize=figsize)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=25)
    if grid:
        ax.grid()

    return fig, ax


BOLSIG_data = np.loadtxt('swarm.txt', skiprows=41, max_rows=100)

# Saving to a .mat file
bolsig_output_arry = BOLSIG_data[:,[1,3,4,11,12]]
mdic = {'BOLSIG_data': bolsig_output_arry, 'label': 'data'}
savemat("bolsig_data.mat", mdic)

BOLSIG_EN = BOLSIG_data[:,1]
BOLSIG_mu = BOLSIG_data[:,3]
BOLSIG_D = BOLSIG_data[:,4]
BOLSIG_alpha = BOLSIG_data[:,11]
BOLSIG_eta = BOLSIG_data[:,12]


# Fit mobility coefficient - electrons
# poly_model = odr.polynomial(1)
# data = odr.Data(np.log10(BOLSIG_EN), np.log10(BOLSIG_mu))
# odr_obj = odr.ODR(data, poly_model)
# output = odr_obj.run()  # running ODR fitting
# poly = np.poly1d(output.beta[::-1])
# print(output.beta[::-1])
# mu_reg = poly(np.log10(BOLSIG_EN))

# plt.loglog(BOLSIG_EN, BOLSIG_mu)
# plt.loglog(BOLSIG_EN, np.power(10, mu_reg), lw=5)
# plt.loglog(BOLSIG_EN, 10**24.76713228 * 10**(-0.34808749*np.log10(BOLSIG_EN)))
# plt.show()
# exit()

# Result: fit coeffs to linear eq are b=24.76713228, m=-0.34808749


# Fit diffusion coefficient - electrons
# poly_model = odr.polynomial(1)

# x = np.log10(BOLSIG_EN[:72])
# y=BOLSIG_D[:72]
# coeffs1 = np.polyfit(x,y,1)
# poly = np.poly1d(coeffs1)
# print(f'First curve: m={coeffs1[0]}, b={coeffs1[1]}')
# diffusion_reg1 = poly(x)

# x = np.log10(BOLSIG_EN[72:])
# y=BOLSIG_D[72:]
# coeffs2 = np.polyfit(x,y,1)
# poly = np.poly1d(coeffs2)
# print(f'Second curve: m={coeffs2[0]}, b={coeffs2[1]}')
# diffusion_reg2 = poly(x)

# log_EN = np.log10(BOLSIG_EN)
# step = 0.5*(np.tanh(1000*(log_EN-1.9))+1)
# diffusion_combined = (1-step)*(3.477702502838188e+23*log_EN + 1.4137092744995724e+24)  + step*(4.6338371114205166e+24*log_EN + -6.700654851204797e+24)

# plt.semilogx(BOLSIG_EN, BOLSIG_D)
# plt.semilogx(BOLSIG_EN[:72], diffusion_reg1, lw=5)
# plt.semilogx(BOLSIG_EN[72:], diffusion_reg2, lw=5)
# plt.semilogx(BOLSIG_EN, diffusion_combined, 'black')
# plt.show()
# exit()


# Fit ionization coefficients
# Ionization
log_EN = np.log10(BOLSIG_EN)

x2 = np.log10(BOLSIG_EN[57:])
y2 = np.log10(BOLSIG_alpha[57:])
coeffs_alpha = np.polyfit(x2,y2,4)
print(f'Ionization curve: a0={coeffs_alpha[4]}, a1={coeffs_alpha[3]}, a2={coeffs_alpha[2]}, a3={coeffs_alpha[1]}, a4={coeffs_alpha[0]}')
poly = np.poly1d(coeffs_alpha)
alpha_reg = poly(log_EN)

step = 0.5*(np.tanh(1000*(BOLSIG_EN-20.09))+1)
# alpha_combined = step*np.power(10, alpha_reg)
alpha_log = -212.3274758440036 + 304.25126942475583*log_EN -186.28749042736393*log_EN**2 + 51.50075493521921*log_EN**3 -5.3618795671332915*log_EN**4
alpha_combined = step*np.power(10, alpha_log)


# plt.loglog(BOLSIG_EN, BOLSIG_alpha)
# plt.loglog(BOLSIG_EN, alpha_combined)
# plt.show()


# Attachment
log_EN = np.log10(BOLSIG_EN)
y_all = np.log10(BOLSIG_eta)

poly_model = odr.polynomial(1)
x1 = np.log10(BOLSIG_EN[:52])
y1 = np.log10(BOLSIG_eta[:52])
coeffs1 = np.polyfit(x1,y1,1)
print(f'First curve: m={coeffs1[0]}, b={coeffs1[1]}')
poly = np.poly1d(coeffs1)
eta_reg1 = poly(log_EN)

x2 = np.log10(BOLSIG_EN[52:])
y2 = np.log10(BOLSIG_eta[52:])
coeffs_eta = np.polyfit(x2,y2,4)
print(f'Second curve: a0={coeffs_eta[4]}, a1={coeffs_eta[3]}, a2={coeffs_eta[2]}, a3={coeffs_eta[1]}, a4={coeffs_eta[0]}')
poly = np.poly1d(coeffs_eta)
eta_reg2 = poly(log_EN)

step = 0.5*(np.tanh(1000*(log_EN-1.1))+1)
# eta_combined = (1-step)*eta_reg1 + step*eta_reg2
eta_reg1 = -40.12938813577186 -1.044327759364851*log_EN
eta_reg2 = -125.18044468718486 + 168.33587943462615*log_EN -103.64966362614207*log_EN**2 + 28.456906756759984*log_EN**3 -2.9427386149316104*log_EN**4
eta_combined = (1-step)*eta_reg1 + step*eta_reg2

plt.loglog(BOLSIG_EN, BOLSIG_eta)
plt.loglog(BOLSIG_EN, np.power(10, eta_combined))
plt.loglog(BOLSIG_EN, BOLSIG_alpha)
plt.loglog(BOLSIG_EN, alpha_combined)
plt.show()


# Next, plot with the raw coefficients here: