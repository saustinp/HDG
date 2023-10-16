import numpy as np
import matplotlib.pyplot as plt

z_dim = np.linspace(0,1e-3,500)
z_nondim = np.linspace(0,2,500)

De = 0.12
t0 = 2e-9
vd = 1.7e5
r_tip = 5e-4
r=0
net_alpha = 5009.5
t0_nondim = t0/(r_tip/vd)

u0_dim = lambda z, t: (4*np.pi*De*t)**(-3/2) * np.exp((z-vd*t)**2/(-4*De*t)  +  net_alpha*vd*t)
u0_nondim = lambda z,t: (4*np.pi*De*t/r_tip**2* (r_tip/vd))**(-3/2) * np.exp(((z-t)**2)/(-4*De*t/r_tip**2* r_tip/vd) + net_alpha*vd*t * r_tip/vd)

max_dim = u0_dim(vd*t0, t0)
max_nondim = u0_nondim(t0_nondim, t0_nondim)
print(max_dim)
print(max_nondim)

print('Check: these should be equal:')
print(max_dim/max_nondim)
print(r_tip**-3)

fig, (ax1, ax2) = plt.subplots(2,1)
ax1.plot(z_dim, u0_dim(z_dim, t0))
ax2.plot(z_nondim, u0_nondim(z_nondim, t0_nondim))
plt.show()