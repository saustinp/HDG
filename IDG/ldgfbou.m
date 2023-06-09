function fh = ldgfbou(ib,ui,nl,p,udg,uh,param,time)
    % fhg = fbou(ib, uinf, nlg1, pg1, udg1, uhg, app.arg, app.time);
    %FHAT flux function
    %   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
    %
    %      NL(N,ND)              Normal N points
    %      P(N,ND)               Coordinates for N points
    %      U(N,NC)               Unknown vector for N points with NC components
    %      Q(N,NC,ND)            Flux vector for N points with NC components in the
    %                            coordinate directions
    %      M(N,NC)               Hybrid unkowns for N points with NC components
    %      PARAM                 Parameter list
    %      FH(N,NC):              Volume flux at N points
    %      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
    %      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
    %      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q
    
    % if ib == 1                 % Far field
    %     um = repmat(ui,size(udg,1),1);
    % elseif  ib == 2            % Inviscid wall
    %     un = udg(:,2).*nl(:,1)+udg(:,3).*nl(:,2);
    %     um = [udg(:,1),udg(:,2)-2*un.*nl(:,1),udg(:,3)-2*un.*nl(:,2),udg(:,4)];
    % end 
    % 
    % %fn = ldgfhat(up,um,np,p,param,time);
    % fh = ldgfhat(nl,p,udg,um,uh,param,time);
    
    up = udg;
    np = nl;
    if ib == 1                 % Far field
        um = repmat(ui,size(up,1),1);
        fh = ldgfhat(nl,p,up,um,uh,param,time);
    elseif  ib == 2            % Reflect
        un = up(:,2).*np(:,1)+up(:,3).*np(:,2);
        um = [up(:,1),up(:,2)-2*un.*np(:,1),up(:,3)-2*un.*np(:,2),up(:,4)];
        fh = ldgfhat(nl,p,up,um,uh,param,time);
    %     uv = up(:,2,:)./up(:,1,:);
    %     vv = up(:,3,:)./up(:,1,:);
    %     gam = param{1};
    %     p = (gam-1)*(up(:,4,:) - 0.5*(up(:,2,:).*uv + up(:,3,:).*vv));
    %     fh = 0*up;
    %     fh(:,2) = p.*nl(:,1);
    %     fh(:,3) = p.*nl(:,2);
    %     [f1 fh]
    %     pause
    end 
    %fh = ldgfhat(nl,p,up,um,uh,param,time);
    
    %fh = euleri_roe(up,um,np,p,param,time);
    
    % %[ng,nc] = size(udg);
    % nch = 1;
    % 
    % tau   = param{end};
    % u     = udg(:,nch);
    % switch ib
    %     case 1  % Dirichlet                
    %         c  = param{2};
    %         qn = udg(:,2).*nl(:,1) + udg(:,3).*nl(:,2);        
    %         fh = qn + (nl*c(:)).*u + tau*(u-ui);                
    %     case 2  
    %         fh = ui;        
    %     case 3  % Prescribed flux
    % %         x = p(:,1);
    % %         y = p(:,2);
    % %         ui = sin(x)*sin(y);            
    %     otherwise
    %         error('unknown boundary type');
    % end
    % 
    % 
    % % [ng,nc] = size(udg);
    % % nch = 1;
    % % nq = nc-nch;
    % % nd = nq;
    % % 
    % % kappa = param{1};
    % % tau   = param{end};
    % % 
    % % u1 = udg1(:,1);
    % % q1x = udg1(:,2);
    % % q1y = udg1(:,3);
    % % u2 = udg2(:,1);
    % % q2x = udg2(:,2);
    % % q2y = udg2(:,3);
    % % qx = 0.5*(q1x+q2x);
    % % qy = 0.5*(q1y+q2y);
    % % 
    % % % ng x nch
    % % fh = kappa.*(qx.*nl(:,1)+qy.*nl(:,2)) + tau.*(u1(:,1)-u2(:,1));
    % % 
    % % 
    % % 
    