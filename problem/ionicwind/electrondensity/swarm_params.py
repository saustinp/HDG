import numpy as np
import matplotlib.pyplot as plt

def get_swarm_params(E, N):
    # Using the values in the appendix of this paper: https://iopscience.iop.org/article/10.1088/0022-3727/30/4/017/pdf
    # Note that the reduced E field is converted to V*cm2 for performing the calculations, then the resultant quantity is converted back from cm->m.

    cm2m = 0.01
    m2cm = 100
    E_N = E/N*m2cm**2   # Reduced electric field in V*cm2. Expects E to be the L2 norm of the E field.

    # First ionization coefficient
    if E_N <= 1.5e-15:
        alpha = 6.619e-17*np.exp(-5.593e-15/E_N)

    elif E_N > 1.5e-15:
        alpha = 2e-16*np.exp(-7.248e-15/E_N)

    alpha *= cm2m**2      # Unit conversion

    # Two-body attachment coefficient
    if E_N <= 1.05e-15:
        eta2 = 6.089e-4*E_N - 2.893e-19

    elif E_N > 1.05e-15:
        eta2 = 8.889e-5*E_N + 2.567e-19

    eta2 *= cm2m**2      # Unit conversion

    # Recombination coefficient
    beta = 2e-7*cm2m**3

    # Electron mobility
    # Removing the sign term because these formulas were designed for computing the velocity, not the mobility. The mobility always has the same sign.
    # Also removed the sign, need to be unsigned for the simulations.

    if E_N <= 2.6e-17:
        mue = (3.38e4 + 6.87e22*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility
    
    elif E_N > 2.6e-17 and E_N <= 1e-16:
        mue = (1.63e6 + 7.2973e21*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility

    elif E_N > 1e-16 and E_N <= 2e-15:
        mue = (1.3e6 + 1.03e22*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility

    elif E_N > 2e-15:
        mue = (7.1e6 + 7.4e21*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility
    
    # Negative ion mobility
    if E_N <= 5e-16:
        mun = 1.86e-4
    elif E_N > 5e-16:
        mun = 2.7e-4

    # 2e-4 m^2/v/s for ion mobility ^


    # Positive ion mobility
    # Assuming that P0/P is approximately 1 and that the units for the E field are V/cm. For now assume that P is approx equal to P0.
    mup = 2.34e-4

    # Electron diffusion coefficient
    # Unit conversion (cm->m) is "baked in" through the mue calculated previously
    D = mue*0.3341e9*E_N**0.54069

    return alpha, eta2, beta, D, mue, mup, mun, E_N/m2cm**2

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

"""
Calc the neutral number density from the ideal gas law
PV=Nk_bT    https://en.wikipedia.org/wiki/Ideal_gas_law


"""

BOLSIG_data = np.loadtxt('swarm.txt', skiprows=41, max_rows=100)
# Need to pull
# A2, A6, A18, A19
BOLSIG_EN = BOLSIG_data[:,1]
BOLSIG_mobility = BOLSIG_data[:,3]
BOLSIG_diffusion = BOLSIG_data[:,4]
BOLSIG_ionization = BOLSIG_data[:,11]
BOLSIG_attach = BOLSIG_data[:,12]

P = 101325 # Pa
V = 1 # m^3
T = 298.15 # K
k_b = 1.380649e-23 #m2 kg s-2 K-1
N = P*V/(k_b*T)
rho = 101325

print(f'N: {N} particles/m^3')

# E_bd = 3e8
# print(f'mue at E_bd: {mue}')
# print(f'D at E_bd: {D}')
# print(f'alpha at E_bd: {alpha}')
# print(f'eta at E_bd: {eta2}')

E_vec = np.logspace(3,8,500)

alpha_vec = np.zeros_like(E_vec)
eta2_vec = np.zeros_like(E_vec)
beta_vec = np.zeros_like(E_vec)
D_vec = np.zeros_like(E_vec)
mue_vec = np.zeros_like(E_vec)
mup_vec = np.zeros_like(E_vec)
mun_vec = np.zeros_like(E_vec)
EN_vec = np.zeros_like(E_vec)

for i, E in enumerate(E_vec):
    alpha, eta2, beta, D, mue, mup, mun, EN = get_swarm_params(E, N)

    alpha_vec[i] = alpha
    eta2_vec[i] = eta2
    beta_vec[i] = beta
    D_vec[i] = D
    mue_vec[i] = mue
    mup_vec[i] = mup
    mun_vec[i] = mun
    EN_vec[i] = EN


# fig, ax = create_default_fig('alpha', 'E/N [Td]', 'alpha')
# ax.loglog(EN_vec/1e-21, alpha_vec, label='alpha')
# ax.loglog(BOLSIG_EN, BOLSIG_ionization, label='BOLSIG+')
# ax.legend()
# plt.ylim([1e-30, 1e-19])
# plt.show()

# fig, ax = create_default_fig('eta', 'E/N [Td]', 'eta')
# # ax.loglog(EN_vec/1e-21, eta2_vec, label='eta')
# ax.loglog(BOLSIG_EN, BOLSIG_ionization, label='BOLSIG+')
# ax.loglog(BOLSIG_EN, BOLSIG_attach, label='BOLSIG+')
# # ax.loglog(EN_vec/1e-21, alpha_vec, label='alpha')
# ax.loglog(BOLSIG_EN, BOLSIG_ionization-BOLSIG_attach, label='alpha-eta')
# ax.legend()
# plt.ylim([1e-43, 1e-19])
# plt.show()


# Choose static values at 130 Td for now (index #77 on the list, 76 indexing from 0)
idx=76
alpha_static = BOLSIG_ionization[idx]*N
eta_static = BOLSIG_attach[idx]*N
diffusion_static = BOLSIG_diffusion[idx]/N
mobility_static = BOLSIG_mobility[idx]/N

print('At 130 Td:')
print(f'alpha static: {alpha_static}')
print(f'eta static: {eta_static}')
print(f'diffusion static: {diffusion_static}')
print(f'mobility static: {mobility_static}')

# fig, ax = create_default_fig('alpha', 'E/N [Td]', 'alpha')
# ax.loglog(BOLSIG_EN, BOLSIG_diffusion, label='BOLSIG+')
# ax.legend()
# plt.show()

# fig, ax = create_default_fig('alpha', 'E/N [Td]', 'alpha')
# ax.loglog(BOLSIG_EN, BOLSIG_mobility, label='BOLSIG+')
# ax.legend()
# plt.show()



# Validate against BOLSIG
# Ask Lee for the attachment/recombination coeffs
