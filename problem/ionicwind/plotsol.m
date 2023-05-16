function t = plotsol(i, shwmsh)
% i is the timestep to plot

% https://www.mathworks.com/matlabcentral/answers/96201-how-do-i-access-a-base-workspace-variable-from-within-a-function
mesh = evalin('base', 'mesh');
UDG_history = evalin('base', 'UDG_history');

ne = UDG_history(:,1,:,i);

clf;
txt = sprintf('ne, timestep %d', i);
t = figure(1);
scaplot(mesh,ne,[0,1.2],0,shwmsh); axis equal; axis tight; colormap jet;
xlim([-2,17]);
ylim([-15,7]);
title(txt, 'FontWeight', 'bold');


% ne = sol(:,1,:,i);
% phi = sol(:,2,:,i);
% Er = sol(:,4,:,i);
% Ez = sol(:,6,:,i);
% Emag = sqrt(Er.^2+Ez.^2);

% clf;
% t = tiledlayout(3,2);
% txt = sprintf('Timestep %d', i);
% title(t,txt, 'FontWeight', 'bold')
% nexttile
% % scaplot(mesh,ne,[min(min(ne)) max(max(ne))],0,shwmsh); title('ne');
% scaplot(mesh,ne,[0 0.01327],0,shwmsh); title('ne');
% nexttile
% scaplot(mesh,phi,[min(min(phi)) max(max(phi))],0,shwmsh); title('Phi');
% nexttile
% % scaplot(mesh,Er,[min(min(Er)) max(max(Er))],0,shwmsh); title('Er');
% scaplot(mesh,Er,[-1000*220e-6, 1000*220e-6],0,shwmsh); title('Er');
% nexttile
% % scaplot(mesh,Ez,[min(min(Ez)) max(max(Ez))],0,shwmsh); title('Ez');
% scaplot(mesh,Ez,[0 2000*220e-6],0,shwmsh); title('Ez');
% nexttile
% ndotE = -Ez;    % n=(0,-1)
% alpha = 0.5*(tanh(1000000*ndotE)+1);
% scaplot(mesh,alpha,[0 1],0,shwmsh); title('alpha');
% % xlim([-.005, 0.03]);
% % ylim([-.055,-0.03]);
% nexttile
% scaplot(mesh,Emag,[0, 1000*220e-6],0,shwmsh); title('E mag');


end