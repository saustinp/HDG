% Plotting the E field
Er0 = mesh.dgnodes(:,3,:);
Ez0 = mesh.dgnodes(:,4,:);
Er_prime = UDG_history(:,8,:,itime+1);
Ez_prime = UDG_history(:,12,:,itime+1);
Er = Er_prime + Er0;
Ez = Ez_prime + Ez0;
normE = sqrt(Er.^2 + Ez.^2);
N = 2.4614924955148245e25;

logEN_Td = log10(normE*E_bd/N/1e-21);

maskEN = normE*E_bd/N/1e-21>130;

scaplot(mesh,maskEN,[-.001, 1.001],0,0); axis equal; axis tight; colormap jet; title('');