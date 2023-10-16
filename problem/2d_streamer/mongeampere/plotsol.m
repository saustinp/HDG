function plotsol(mesh, master, UDG, ACG, str, href)

Mach = 7;
mf1 = 18;

figure(1); clf; 
meshplot(mesh,1);
set(gca,'FontSize',mf1); 
set(gca,'LooseInset',get(gca,'TightInset'));
axis on; box on; axis tight;
ax = gca;
fn = str + "_mesh" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

figure(2); clf; 
scaplot(mesh, ACG,[],2); colormap jet;
set(gca,'FontSize',mf1); 
axis off; box on; axis tight;
ax = gca;
fn = str + "_solution_av" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

figure(3); clf; 
scaplot(mesh, eulereval(UDG,'M',1.4,Mach),[0 7],2); colormap jet;
set(gca,'FontSize',mf1); 
axis off; box on; axis tight;
ax = gca;
fn = str + "_solution_Mach" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

figure(4); clf; 
scaplot(mesh, eulereval(UDG,'r',1.4,Mach),[1 5.7],2); colormap jet;
set(gca,'FontSize',mf1); 
axis off; box on; axis tight;
ax = gca;
fn = str + "_solution_density" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

figure(5); clf; 
scaplot(mesh, eulereval(UDG,'p',1.4,Mach),[0.015 0.92],2); colormap jet;
set(gca,'FontSize',mf1); 
axis off; box on; axis tight;
ax = gca;
fn = str + "_solution_pressure" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

if nargin>5
s = sensor(mesh, master, UDG, href, mesh.dist);
figure(6); clf; 
scaplot(mesh, eulereval(s,'r',1.4,Mach),[],2); colormap jet;
set(gca,'FontSize',mf1); 
axis off; box on; axis tight;
ax = gca;
fn = str + "_solution_sensor" + ".png";
exportgraphics(ax,fn,'Resolution',200); 
end

