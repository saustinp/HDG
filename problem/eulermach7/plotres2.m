% clear p1 r1;
% for i = 1:length(UDG1)
%     [x1,y1,p1(:,i)] = getfieldatbou(mesh,eulereval(UDG1{i},'p',1.4,7));
%     [x1,y1,r1(:,i)] = getfieldatbou(mesh,eulereval(UDG1{i},'r',1.4,7));
% end
% ii = [1 3 5 7];
% figure(1); clf; hold on;
% plot(y1,p1(:,ii),'-','LineWidth',2,'MarkerSize',8);
% figure(2); clf; hold on;
% plot(y1,r1(:,ii),'-','LineWidth',2,'MarkerSize',8);

% clear p2 r2;
% for i = 1:length(UDG2)
%     [x2,y2,p2(:,i)] = getfieldatbou(mesh2,eulereval(UDG2{i},'p',1.4,7));
%     [x2,y2,r2(:,i)] = getfieldatbou(mesh2,eulereval(UDG2{i},'r',1.4,7));
% end

for i = 1:12
  figure(1); clf;
  scaplot(mesh2,eulereval(UDG2{i},'p',1.4,7),[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; colorbar off;
  ax = gca;
  fn = "cyl7_pressure" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
  figure(2); clf;
  scaplot(mesh2,eulereval(UDG2{i},'r',1.4,7),[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; colorbar off;
  ax = gca;
  fn = "cyl7_density" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
  figure(3); clf;
  scaplot(mesh2,eulereval(UDG2{i},'M',1.4,7),[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; %colorbar off;
  ax = gca;
  fn = "cyl7_mach" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
  figure(4); clf;
  scaplot(mesh2,ACG2{i},[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; colorbar off;
  ax = gca;
  fn = "cyl7_av" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
end

figure(3); clf;
scaplot(mesh2,eulereval(UDG2{i},'H',1.4,7),[],2,0); 
set(gca,'FontSize',16); 
%set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off; %colorbar off;
ax = gca;
fn = "cyl7_enthalpy" + num2str(i) + ".png";
exportgraphics(ax,fn,'Resolution',300); 

  figure(3); clf;
  scaplot(mesh,eulereval(UDG1{i},'H',1.4,7),[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; %colorbar off;
  ax = gca;
  fn = "cyl7_coarseenthalpy" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 


 figure(3); clf;
  scaplot(mesh,eulereval(UDG1{i},'M',1.4,7),[0 3],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; %colorbar off;
  ax = gca;
  fn = "cyl3_coarsemach" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
figure(3); clf;
scaplot(mesh2,eulereval(UDG2{i},'M',1.4,7),[],2,0); 
set(gca,'FontSize',16); 
%set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off; %colorbar off;
ax = gca;
fn = "cyl3_mach" + num2str(i) + ".png";
exportgraphics(ax,fn,'Resolution',300); 
  
for i = 1:12    
  figure(4); clf;
  scaplot(mesh2,ACG2{i},[],2,0); 
  set(gca,'FontSize',16); 
  axis equal; axis tight; axis off;   
  fn = "cyl7_avcolorbar" + num2str(i) + ".png";  
  ax = gca;
  exportgraphics(ax,fn,'Resolution',300); 
end

for i = 1:9
%   figure(1); clf;
%   scaplot(mesh,eulereval(UDG1{i},'p',1.4,7),[],2,0); 
%   set(gca,'FontSize',16); 
%   %set(gca,'LooseInset',get(gca,'TightInset'))
%   axis equal; axis tight; axis off; colorbar off;
%   ax = gca;
%   fn = "cyl7_coarsepressure" + num2str(i) + ".png";
%   exportgraphics(ax,fn,'Resolution',300); 
  
  figure(2); clf;
  scaplot(mesh,eulereval(UDG1{i},'r',1.4,7),[0.9 5.86],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; colorbar off;
  ax = gca;
  fn = "cyl7_coarsedensity" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
  figure(3); clf;
  scaplot(mesh,eulereval(UDG1{i},'M',1.4,7),[0 7],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; %colorbar off;
  ax = gca;
  fn = "cyl7_coarsemach" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
  
  figure(4); clf;
  scaplot(mesh,ACG1{i},[],2,0); 
  set(gca,'FontSize',16); 
  %set(gca,'LooseInset',get(gca,'TightInset'))
  axis equal; axis tight; axis off; %colorbar off;
  ax = gca;
  fn = "cyl7_coarseav" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',300); 
end


ii = [1 5 9 12];

figure(1); clf; hold on;
plot(y2,p2(:,ii),'-','LineWidth',2,'MarkerSize',8);
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight; axis on; box on;
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$\mbox{Pressure}$", 'interpreter', 'latex', 'FontSize', 26);
% leg = legend({"$n=1$","$n=5$","$n=9$","$n=12$"}, 'interpreter', 'latex', 'FontSize', 19, 'Location', 'NE');
% leg.ItemTokenSize = [30,10];
axis([-1 1 0 1]);
%print -dpng cyl7_surfacepressure.png

figure(2); clf; hold on;
plot(y2,r2(:,ii),'-','LineWidth',2,'MarkerSize',8);
set(gca,'FontSize',20); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight; axis on; box on;
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 24);
ylabel("$\mbox{Density}$", 'interpreter', 'latex', 'FontSize', 24);
leg = legend({"$n=1$","$n=5$","$n=9$","$n=12$"}, 'interpreter', 'latex', 'FontSize', 19, 'Location', 'NE');
leg.ItemTokenSize = [30,10];
axis([-1 1 1 6]);
%print -dpng cyl7_surfacedensity.png

figure(3); clf;
scaplot(mesh2,ACG2{1}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.002','0.004','0.006','0.008','0.010','0.012','0.014'}, ...
   'XTick', 0:0.002:0.014);
print -dpng cyl7_av1.png

figure(3); clf;
scaplot(mesh2,ACG2{5}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.001','0.002','0.003','0.004','0.005','0.006'}, ...
   'XTick', 0:0.001:0.006);
print -dpng cyl7_av5.png

figure(3); clf;
scaplot(mesh2,ACG2{7}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.0005','0.0010','0.0015','0.0020','0.0025','0.0030','0.0035','0.0040'}, ...
   'XTick', 0:0.0005:0.004);
print -dpng cyl7_av7.png

figure(3); clf;
scaplot(mesh2,ACG2{9}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.0005','0.001','0.0015','0.0020','0.0025'}, ...
   'XTick', 0:0.0005:0.0025);
print -dpng cyl7_av9.png

figure(3); clf;
scaplot(mesh2,ACG2{12}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.0002','0.0004','0.0006','0.0008','0.0010','0.0012'}, ...
   'XTick', 0:0.0002:0.0012);
print -dpng cyl7_av12.png

figure(3); clf;
scaplot(mesh2,eulereval(UDG2{12},'p',1.4,7),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off; colorbar off;
print -dpng cyl7_pressure12.png


clear rho2 pres2 mach2
for i = 1:length(UDG2)
    [x2,rho2(:,i)] = getfieldaty(mesh2,eulereval(UDG2{i},'r',1.4,7));    
    [x2,pres2(:,i)] = getfieldaty(mesh2,eulereval(UDG2{i},'p',1.4,7));    
    [x2,mach2(:,i)] = getfieldaty(mesh2,eulereval(UDG2{i},'M',1.4,7));    
end

poly = polyfit(pn3(:,2) ,pn3(:,1), 6);
y1 = linspace(-2.645,2.645, 1000);
x1 = polyval(poly,y1);

msh = mesh;
msh.p(:,1) = mesh.p(:,2);
msh.p(:,2) = -mesh.p(:,1);
msh.dgnodes(:,1,:) = mesh.dgnodes(:,2,:);
msh.dgnodes(:,2,:) = -mesh.dgnodes(:,1,:);

figure(1); clf; hold on;
meshplot(mesh,1);
plot(x1,y1,'-r','LineWidth',1.5);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
print -dpng cyl7_mesh1.png


ii=[4 7 10 12];
figure(2); clf; 
plot(x2, rho2(:,ii), '-', 'LineWidth', 2); 
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight; axis on; box on;
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 24);
leg = legend({"$n=4$" , ...
              "$n=7$" , ...
              "$n=10$" , ...
              "$n=12$" , ...              
              }, 'interpreter', 'latex', 'FontSize', 18, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl7_centerlinedensity.png

figure(2); clf; 
plot(x2, pres2(:,ii), '-', 'LineWidth', 2); 
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight; axis on; box on;
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
ylabel("$p$", 'interpreter', 'latex', 'FontSize', 24);
leg = legend({"$n=4$" , ...
              "$n=7$" , ...
              "$n=10$" , ...
              "$n=12$" , ...              
              }, 'interpreter', 'latex', 'FontSize', 18, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl7_centerlinepressure.png

ii = [4 7 10 12];
figure(2); clf; 
plot(x2, mach2(:,ii), '-', 'LineWidth', 2); 
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight; axis on; box on;
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
ylabel("$M$", 'interpreter', 'latex', 'FontSize', 24);
leg = legend({"$n=4$" , ...
              "$n=7$" , ...
              "$n=10$" , ...
              "$n=12$" , ...              
              }, 'interpreter', 'latex', 'FontSize', 18, 'Location', 'NE');
leg.ItemTokenSize = [25,10];
print -dpng cyl7_centerlinemach.png
