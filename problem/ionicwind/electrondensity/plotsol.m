function t = plotsol(i, shwmsh)
    % i is the timestep to plot
    
    % https://www.mathworks.com/matlabcentral/answers/96201-how-do-i-access-a-base-workspace-variable-from-within-a-function
    mesh = evalin('base', 'mesh');
    UDG_history = evalin('base', 'UDG_history');
    
    ne = UDG_history(:,1,:,i);
    nn = UDG_history(:,2,:,i);
    np = UDG_history(:,3,:,i);
    phi = UDG_history(:,4,:,i);
    Er = UDG_history(:,8,:,i);
    Ez = UDG_history(:,12,:,i);
    
    clf;
    t = tiledlayout(3,2);
    txt = sprintf('Timestep %d', i);
    title(t,txt, 'FontWeight', 'bold')
    nexttile
    % scaplot(mesh,ne,[min(min(ne)) max(max(ne))],0,shwmsh); title('ne');
    scaplot(mesh,ne,[min(min(ne)) max(max(ne))],0,shwmsh); title('ne');
    xlim([-2,8]);
    ylim([-10,2]);
    nexttile
    scaplot(mesh,nn,[min(min(nn)) max(max(nn))],0,shwmsh); title('nn');
    xlim([-2,8]);
    ylim([-10,2]);
    nexttile
    scaplot(mesh,np,[min(min(np)) max(max(np))],0,shwmsh); title('np');
    xlim([-2,8]);
    ylim([-10,2]);
    nexttile
    scaplot(mesh,phi,[-15.15 0],0,shwmsh); title('Phi');
    xlim([-2,8]);
    ylim([-10,2]);
    nexttile
    % scaplot(mesh,Er,[min(min(Er)) max(max(Er))],0,shwmsh); title('Er');
    scaplot(mesh,Er,[min(min(Er)) max(max(Er))],0,shwmsh); title('Er');
    xlim([-2,8]);
    ylim([-10,2]);
    nexttile
    % scaplot(mesh,Ez,[min(min(Ez)) max(max(Ez))],0,shwmsh); title('Ez');
    scaplot(mesh,Ez,[min(min(Ez)) max(max(Ez))],0,shwmsh); title('Ez');
    xlim([-2,8]);
    ylim([-10,2]);
    
    fname = sprintf('./out_frames/barchart%d.png', i);
    exportgraphics(t,fname,'Resolution',300);

    
    end