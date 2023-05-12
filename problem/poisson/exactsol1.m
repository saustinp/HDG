function uex = exactsol1(p)

x = p(:,1,:);
y = p(:,2,:);
uex(:,1,:) = sin(0.5*pi*x).*sin(0.5*pi*y);       
uex(:,2,:) = -(pi*cos((pi*x)/2).*sin((pi*y)/2))/2;
uex(:,3,:) = -(pi*cos((pi*y)/2).*sin((pi*x)/2))/2;
 
