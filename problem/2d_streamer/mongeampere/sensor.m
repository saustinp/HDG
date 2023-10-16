function s = sensor(mesh, master, UDG, href, dist, beta)

if nargin<6
  beta = 1;
end

  
% r    = UDG(:,1,:);
% ru   = UDG(:,2,:);
% rv   = UDG(:,3,:);
% rE   = UDG(:,4,:);
% rx   = UDG(:,5,:);
% rux  = UDG(:,6,:);
% rvx  = UDG(:,7,:);
% rEx  = UDG(:,8,:);
% ry   = UDG(:,9,:);
% ruy  = UDG(:,10,:);
% rvy  = UDG(:,11,:);
% rEy  = UDG(:,12,:);

% gam = 1.4;
% gam1 = gam-1;
% r1   = 1./r;
% uv   = ru.*r1;
% vv   = rv.*r1;
% q   = 0.5*(uv.*uv+vv.*vv);
% p    = gam1*(rE-r.*q);
                                        
% ux  = (rux - rx.*uv).*r1;
% vx  = (rvx - rx.*vv).*r1;
% qx  = uv.*ux + vv.*vx;
% px  = gam1*(rEx - rx.*q - r.*qx);

% uy  = (ruy - ry.*uv).*r1;
% vy  = (rvy - ry.*vv).*r1;
% qy  = uv.*uy + vv.*vy;
% py  = gam1*(rEy - ry.*q - r.*qy);
 
% alpha = 1000;    
% s = (rx.^2 + ry.^2);
% s = limiting(s,0,max(s(:))/2,alpha,0);
% if isempty(dist)
%   s = sqrt(beta + s);
% else
%   s = sqrt(beta + s.*dist);  
% end

Er = UDG(:,6,:);
Ez = UDG(:,9,:);
normE = sqrt(Er.^2+Ez.^2);  % This will be used as the sensor field. Note not multiplying by 3e6 because it will just be normalized anyway on the next line

% normalization
normE = normE/max(max(normE))*10;   % Rescaling to [0,10]
s = sqrt(beta+normE);
  
%min(s(:))

mesh.ib = [];
mesh.in = 1:size(mesh.p2,1);
s = cgpoisson(mesh, master, s, [href 1.0]);    
s = s(mesh.t2');
s = reshape(s, size(UDG(:,1,:)));
s = s/min(s(:));
%s = 1 + c*s;

% figure(1); clf; 
% scaplot(mesh, s(:,1,:),[],2,1); axis on; axis equal; axis tight;

% figure(2); clf; 
% scaplot(mesh, py(:,1,:),[],2,1); axis on; axis equal; axis tight;
