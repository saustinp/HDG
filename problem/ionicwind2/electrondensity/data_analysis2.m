itime=1;

param = init_phys_param();
mue_ref = param{19};
r_tip = param{16};
Kep = param{5};
E_bd = param{15};
e = param{11};
epsilon0 = param{12};
    
Er0 = mesh.dgnodes(:,3,:);
Ez0 = mesh.dgnodes(:,4,:);

ne = UDG_history(:,1,:,itime);
nn = UDG_history(:,2,:,itime);
np = UDG_history(:,3,:,itime);
Er_prime = UDG_history(:,8,:,itime);
Ez_prime = UDG_history(:,12,:,itime);

Er = Er_prime + Er0;
Ez = Ez_prime + Ez0;
normE = sqrt(Er.^2 + Ez.^2);
N = 2.4614924955148245e25;

EN_Td = normE*E_bd/N/1e-21;
logEN_Td = log10(EN_Td);

EN_Td_thresh = normE*E_bd/N/1e-21>700;
EN_Td_thresh = EN_Td_thresh+1;

E_Vm = normE*E_bd;

% er = UDG_history(:,8,:,1)*E_bd;
% ez = UDG_history(:,12,:,1)*E_bd;
% normE = sqrt(er.^2 + ez.^2);
% EN_Td = normE/N/1e-21;

% figure();
% scaplot(mesh,E_Vm,[],0,0); axis equal; axis tight; colormap jet; title('E
% field magnitude [V/m]'); xlim([-.2, 1.2]); ylim([-.8,.4]);4
max(max(logEN_Td))

figure();
scaplot(mesh,EN_Td,[],0,0); axis equal; axis tight; colormap jet; title('log10(reduced E field magnitude in Td)');% xlim([-.2, 1.2]); ylim([-.8,.4]);


% figure();
% scaplot(mesh,EN_Td_thresh,[.99, 2.001],0,0); axis equal; axis tight; colormap jet; title('Reduced E field in Td > 130'); xlim([-.2, 1.2]); ylim([-.8,.4]);
return;

% Plot alpha
% Should be:
% At farfield (EN_Td = 2, E=5.3736e+04 V/m): ~0
% At farfield (EN_Td = 1000, E=2.6868e+07 V/m): 2.9942e+05

alpha = get_alpha(E_Vm, N);
% figure();
% scaplot(mesh,alpha,[],0,0); axis equal; axis tight; colormap jet; title('alpha'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot eta
% Should be:
% At farfield (EN_Td = 2): 9.6709e-16
% At farfield (EN_Td = 1000): 242.1575

eta = get_eta(E_Vm, N);
% figure();
% scaplot(mesh,eta,[],0,0); axis equal; axis tight; colormap jet; title('eta'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot diffusion
% Should be:
% At farfield (EN_Td = 2): 0.0565
% At farfield (EN_Td = 1000): 0.2680

D = get_diffusion_e(E_Vm, N);
% figure();
% scaplot(mesh,D,[],0,0); axis equal; axis tight; colormap jet; title('diffusion coefficient'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot mu
% Should be:
% At farfield (EN_Td = 2): 0.1710
% At farfield (EN_Td = 1000): 5.2e23/N = 0.0197

mue = get_mue(E_Vm, N);
figure();
scaplot(mesh,mue,[],0,0); axis equal; axis tight; colormap jet; title('mue'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot ne
figure();
scaplot(mesh,ne-np+nn,[],0,0); axis equal; axis tight; colormap jet; title('total'); xlim([-.2, 1.2]); ylim([-.8,.4]);

figure();
scaplot(mesh,np,[],0,0); axis equal; axis tight; colormap jet; title('np'); xlim([-.2, 1.2]); ylim([-.8,.4]);

figure();
scaplot(mesh,ne,[],0,0); axis equal; axis tight; colormap jet; title('ne'); xlim([-.2, 1.2]); ylim([-.8,.4]);

figure();
scaplot(mesh,nn,[],0,0); axis equal; axis tight; colormap jet; title('nn'); xlim([-.2, 1.2]); ylim([-.8,.4]);

figure();
scaplot(mesh,Er_prime,[],0,0); axis equal; axis tight; colormap jet; title('Er_prime'); xlim([-.2, 1.2]); ylim([-.8,.4]);

figure();
scaplot(mesh,Ez_prime,[],0,0); axis equal; axis tight; colormap jet; title('Ez_prime'); xlim([-.2, 1.2]); ylim([-.8,.4]);
return;

% Plot E
figure();
scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet; title('normE'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot each term in the source expression for ne
figure();
se1 = (alpha-eta).*(mue/mue_ref).*r_tip.*ne.*normE;
scaplot(mesh,se1,[],0,0); axis equal; axis tight; colormap jet; title('(alpha-eta)*(mue/mue_ref)*r_{tip}*ne*normE'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% figure();
% se2 = -Kep.*epsilon0./(e.*mue_ref).*ne.*np;
% scaplot(mesh,se2,[],0,0); axis equal; axis tight; colormap jet; title('-Kep*epsilon0/(e*mue_{ref})*ne*np'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% figure();
% se = se1+se2;
% scaplot(mesh,se,[],0,0); axis equal; axis tight; colormap jet; title('Total source'); xlim([-.2, 1.2]); ylim([-.8,.4]);

% Plot convective velocities
cr_e = -(mue/mue_ref).*Er;
cz_e = -(mue/mue_ref).*Ez;
normC = sqrt(cr_e.^2 + cz_e.^2);

% figure();
% scaplot(mesh,cr_e,[],0,0); axis equal; axis tight; colormap jet; title('cr_e'); xlim([-.2, 1.2]); ylim([-.8,.4]);
% 
% figure();
% scaplot(mesh,cz_e,[],0,0); axis equal; axis tight; colormap jet; title('cz_e'); xlim([-.2, 1.2]); ylim([-.8,.4]);
% 
% figure();
% scaplot(mesh,normC,[],0,0); axis equal; axis tight; colormap jet; title('normc'); xlim([-.2, 1.2]); ylim([-.8,.4]);