import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import odr

def create_default_fig(title, xlabel, ylabel, figsize=(8,6), grid=1):
    """""
    Creates a default figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_title(title, fontsize=17)
    if grid:
        ax.grid()

    return fig, ax

# Load the data from the BOLSIG output file
BOLSIG_data = np.loadtxt('swarm.txt', skiprows=41, max_rows=100)

# Selects the output rows that we want - corresponding to the swarm parameters of interest. The selected columns are written to at matlab .mat file for reading into matlab
bolsig_output_arry = BOLSIG_data[:,[1,3,4,11,12]]
mdic = {'BOLSIG_data': bolsig_output_arry, 'label': 'data'}
savemat("bolsig_data.mat", mdic)

BOLSIG_EN = BOLSIG_data[:,1]
BOLSIG_mu = BOLSIG_data[:,3]
BOLSIG_D = BOLSIG_data[:,4]
BOLSIG_alpha = BOLSIG_data[:,11]
BOLSIG_eta = BOLSIG_data[:,12]

############### MAIN SWITCH ###############
flag = 'fit'
# flag = 'validate'

if flag == 'fit':
    # Fit mobility coefficient - electrons
    poly_model = odr.polynomial(1)
    data = odr.Data(np.log10(BOLSIG_EN), np.log10(BOLSIG_mu))
    odr_obj = odr.ODR(data, poly_model)
    output = odr_obj.run()  # running ODR fitting
    poly = np.poly1d(output.beta[::-1])
    print(output.beta[::-1])
    mu_reg = poly(np.log10(BOLSIG_EN))

    fig, ax = create_default_fig('Mobility coefficient curve fit: regression output', 'Reduced E field [Td]', 'Ionization/attachment coefficient')
    ax.loglog(BOLSIG_EN, BOLSIG_mu, label=r'$\mu$: BOLSIG')
    ax.loglog(BOLSIG_EN, np.power(10, mu_reg), lw=2, label=r'$\mu$: curve fit from scipy')
    ax.legend()
    plt.show()

    # Fit diffusion coefficient - electrons
    poly_model = odr.polynomial(1)

    x = np.log10(BOLSIG_EN[:72])
    y=BOLSIG_D[:72]
    coeffs1 = np.polyfit(x,y,1)
    poly = np.poly1d(coeffs1)
    print(f'First curve: m={coeffs1[0]}, b={coeffs1[1]}')
    diffusion_reg1 = poly(x)

    x = np.log10(BOLSIG_EN[72:])
    y=BOLSIG_D[72:]
    coeffs2 = np.polyfit(x,y,1)
    poly = np.poly1d(coeffs2)
    print(f'Second curve: m={coeffs2[0]}, b={coeffs2[1]}')
    diffusion_reg2 = poly(x)

    log_EN = np.log10(BOLSIG_EN)
    step = 0.5*(np.tanh(1000*(log_EN-1.9))+1)
    diffusion_combined = (1-step)*(3.477702502838188e+23*log_EN + 1.4137092744995724e+24)  + step*(4.6338371114205166e+24*log_EN + -6.700654851204797e+24)

    fig, ax = create_default_fig('Diffusion coefficient fit: regression output', 'Reduced E field [Td]', 'Diffusion coefficient')
    ax.semilogx(BOLSIG_EN, BOLSIG_D, label='D: BOLSIG')
    ax.semilogx(BOLSIG_EN, diffusion_combined, label='D: curve fit from scipy')
    ax.legend()
    plt.show()

elif flag == 'validate':

    # Mobility coefficient
    fig, ax = create_default_fig('Mobility coefficient fit: direct coefficient', 'Reduced E field [Td]', 'Mobility coefficient')
    ax.loglog(BOLSIG_EN, BOLSIG_mu, label=r'$\mu$: BOLSIG')
    ax.loglog(BOLSIG_EN, 10**24.76713228 * 10**(-0.34808749*np.log10(BOLSIG_EN)), label=r'$\mu$: curve fit direct coefficients')
    plt.show()

    # Diffusion coefficient
    log_EN = np.log10(BOLSIG_EN)
    step = 0.5*(np.tanh(1000*(log_EN-1.9))+1)
    diffusion_combined = (1-step)*(3.477702502838188e+23*log_EN + 1.4137092744995724e+24)  + step*(4.6338371114205166e+24*log_EN + -6.700654851204797e+24)

    fig, ax = create_default_fig('Diffusion coefficient fit: direct coefficient', 'Reduced E field [Td]', 'Diffusion coefficient')
    ax.semilogx(BOLSIG_EN, BOLSIG_D, label='D: BOLSIG')
    ax.semilogx(BOLSIG_EN, diffusion_combined, label='D: curve fit direct coefficients')
    plt.show()

    # Ionization coefficients
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

    fig, ax = create_default_fig('BOLSIG data, direct curve fits', 'Reduced E field [Td]', 'Ionization/attachment coefficient')
    ax.loglog(BOLSIG_EN, BOLSIG_eta, label=r'$\eta$: BOLSIG')
    ax.loglog(BOLSIG_EN, np.power(10, eta_combined), label=r'$\eta$: curve fit direct coefficients')
    ax.loglog(BOLSIG_EN, BOLSIG_alpha, label=r'$\alpha$: BOLSIG')
    ax.loglog(BOLSIG_EN, alpha_combined, label=r'$\alpha$: curve fit direct coefficients')
    ax.legend()
    plt.show()