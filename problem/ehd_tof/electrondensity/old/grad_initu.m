syms pi De t z vd r_tip net_alpha r n_ref E_bd

% First, write the nondimensional quantities here:



u0 = (4.*pi.*De.*t./r_tip.^2.* (r_tip./vd)).^(-3./2) .* exp(((z-t).^2 + r.^2)./(-4.*De.*t./r_tip.^2.* r_tip./vd) + net_alpha.*vd.*t .* r_tip./vd);

du0_dr = -diff(u0, r)   % Need the negative for q. The result is then printed out to the terminal and copy/pasted into the setup script
du0_dz = -diff(u0, z)


% du0_dr = (r*r_tip*vd*exp(net_alpha*r_tip*t - (r_tip*vd*((t - z)^2 + r^2))/(4*De*t)))/(2*De*t*((4*De*pi*t)/(r_tip*vd))^(3/2))
% du0_dz = -(r_tip*vd*exp(net_alpha*r_tip*t - (r_tip*vd*((t - z)^2 + r^2))/(4*De*t))*(2*t - 2*z))/(4*De*t*((4*De*pi*t)/(r_tip*vd))^(3/2))