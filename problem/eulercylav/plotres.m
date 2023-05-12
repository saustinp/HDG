%load results.mat

alpha = 100; beta = 0.5;
s=linspace(0,1.5,1000);
a = (s-beta).*(atan(alpha*(s-beta))/pi + 0.5) - atan(alpha)/pi + 0.5;
figure(1);clf;
plot(s,a,'-','LineWidth',1.5); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
xlabel("$\mbox{Shock strength } S$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\mbox{Source term } s$", 'interpreter', 'latex', 'FontSize', 22);


y1 = linspace(pn(1,2),pn(end,2), 1000);
x1 = polyval(poly,y1);

msh = mesh;
msh.p(:,1) = mesh.p(:,2);
msh.p(:,2) = -mesh.p(:,1);
msh.dgnodes(:,1,:) = mesh.dgnodes(:,2,:);
msh.dgnodes(:,2,:) = -mesh.dgnodes(:,1,:);

figure(1); clf; hold on;
meshplot(msh,1);
plot(y1,-x1,'-r','LineWidth',1.5);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
print -dpng cyl3_mesh1.png

msh2 = mesh2;
msh2.p(:,1) = mesh2.p(:,2);
msh2.p(:,2) = -mesh2.p(:,1);
msh2.dgnodes(:,1,:) = mesh2.dgnodes(:,2,:);
msh2.dgnodes(:,2,:) = -mesh2.dgnodes(:,1,:);

figure(2); clf; hold on;
meshplot(msh2,1);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
print -dpng cyl3_mesh2.png
 
figure(3); clf;
scaplot(msh,ACG1{3}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.002','0.004','0.006','0.008','0.001'}, ...
   'XTick', 0:0.002:0.01, 'Location', 'northoutside');
print -dpng cyl3_av1.png

figure(3); clf;
scaplot(msh2,ACG2{7}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('XTickLabel',{'0','0.0001','0.0002','0.0003','0.0004','0.0005','0.0006','0.0007'}, ...
   'XTick', 0:0.0001:0.0007, 'Location', 'northoutside');
print -dpng cyl3_av2.png

figure(3); clf;
scaplot(msh,UDG1{3}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_density1.png

figure(3); clf;
scaplot(msh2,UDG2{7}(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_density2.png

figure(3); clf;
scaplot(msh,eulereval(UDG1{3},'M',1.4,3),[0 3.025],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_mach1.png

figure(3); clf;
scaplot(msh2,eulereval(UDG2{7},'M',1.4,3),[0 3.025],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_mach2.png

figure(3); clf;
scaplot(msh,eulereval(UDG1{3},'p',1.4,3),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_pres1.png

figure(3); clf;
scaplot(msh2,eulereval(UDG2{7},'p',1.4,3),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis off;
colorbar off;
colorbar('Location', 'northoutside');
print -dpng cyl3_pres2.png

x = linspace(-2.0, -1.0, 1200);
y = [0];
xy = [x(:) y*ones(length(x),1)];
rho1 = fieldatx(mesh,UDG1{3}(:,1,:),xy,7);    
pres1 = fieldatx(mesh,eulereval(UDG1{3},'p',1.4,3),xy,7);    
mach1 = fieldatx(mesh,eulereval(UDG1{3},'M',1.4,3),xy,7);    
rho2 = fieldatx(mesh2,UDG2{7}(:,1,:),xy,10);    
pres2 = fieldatx(mesh2,eulereval(UDG2{7},'p',1.4,3),xy,10);    
mach2 = fieldatx(mesh2,eulereval(UDG2{7},'M',1.4,3),xy,10);    

figure(6); clf; hold on;
plot(-x, mach1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, mach2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$M$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_mach3.png

figure(6); clf; hold on;
plot(-x, pres1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, pres2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$P$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_pres3.png

figure(6); clf; hold on;
plot(-x, rho1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, rho2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_density3.png

x = linspace(-2.0, -0.866, 1200);
y = [0.5];
xy = [x(:) y*ones(length(x),1)];
rho1 = fieldatx(mesh,UDG1{3}(:,1,:),xy,7);    
pres1 = fieldatx(mesh,eulereval(UDG1{3},'p',1.4,3),xy,7);    
mach1 = fieldatx(mesh,eulereval(UDG1{3},'M',1.4,3),xy,7);    
rho2 = fieldatx(mesh2,UDG2{7}(:,1,:),xy,10);    
pres2 = fieldatx(mesh2,eulereval(UDG2{7},'p',1.4,3),xy,10);    
mach2 = fieldatx(mesh2,eulereval(UDG2{7},'M',1.4,3),xy,10);    

figure(6); clf; hold on;
plot(-x, mach1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, mach2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$M$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_mach4.png

figure(6); clf; hold on;
plot(-x, pres1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, pres2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$P$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_pres4.png

figure(6); clf; hold on;
plot(-x, rho1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, rho2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_density4.png

x = linspace(-2.0, 0, 1200);
y = [1.0];
xy = [x(:) y*ones(length(x),1)];
rho1 = fieldatx(mesh,UDG1{3}(:,1,:),xy,7);    
pres1 = fieldatx(mesh,eulereval(UDG1{3},'p',1.4,3),xy,7);    
mach1 = fieldatx(mesh,eulereval(UDG1{3},'M',1.4,3),xy,7);    
rho2 = fieldatx(mesh2,UDG2{7}(:,1,:),xy,10);    
pres2 = fieldatx(mesh2,eulereval(UDG2{7},'p',1.4,3),xy,10);    
mach2 = fieldatx(mesh2,eulereval(UDG2{7},'M',1.4,3),xy,10);    

figure(6); clf; hold on;
plot(-x, mach1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, mach2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$M$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_mach5.png

figure(6); clf; hold on;
plot(-x, pres1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, pres2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$P$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_pres5.png

figure(6); clf; hold on;
plot(-x, rho1, '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(-x, rho2, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng cyl3_density5.png

x = linspace(-2.0, -1, 1200);
y = [0.0];
xy = [x(:) y*ones(length(x),1)];
for i = 1:length(UDG1)
    rho1(:,i) = fieldatx(mesh,UDG1{i}(:,1,:),xy,7);    
    pres1(:,i) = fieldatx(mesh,eulereval(UDG1{i},'p',1.4,3),xy,7);    
    mach1(:,i) = fieldatx(mesh,eulereval(UDG1{i},'M',1.4,3),xy,7);    
end
for i = 1:length(UDG2)
    rho2(:,i) = fieldatx(mesh2,UDG2{i}(:,1,:),xy,10);    
    pres2(:,i) = fieldatx(mesh2,eulereval(UDG2{i},'p',1.4,3),xy,10);    
    mach2(:,i) = fieldatx(mesh2,eulereval(UDG2{i},'M',1.4,3),xy,10);    
end

av = [0.0250    0.0125    0.0063    0.0031    0.0016    0.0008    0.0004];

figure(6); clf; hold on;
plot(-x, mach2(:,2:8), '-', 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$M$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'$\lambda = 0.025$','$\lambda = 0.025/2$', '$\lambda = 0.025/4$', '$\lambda = 0.025/8$','$\lambda = 0.025/16$','$\lambda = 0.025/32$', '$\lambda = 0.025/64$'}, 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [30,10];
print -dpng cyl3_mach6.png

figure(6); clf; hold on;
plot(-x, pres2(:,2:8), '-', 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$P$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'$\lambda = 0.025$','$\lambda = 0.025/2$', '$\lambda = 0.025/4$', '$\lambda = 0.025/8$','$\lambda = 0.025/16$','$\lambda = 0.025/32$', '$\lambda = 0.025/64$'}, 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [30,10];
print -dpng cyl3_density6.png


% figure(4); clf;
% scaplot(mesh,UDG1{3}(:,1,:),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% % xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% % ylabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_sol1.png

% yu = linspace(0.04326, 0.7172, 101);
% xu = polyval(poly,yu);
% 
% yl = linspace(-0.16794,-0.0597, 101);
% xl = polyval(poly2,yl);
% 
% figure(1); clf; hold on;
% meshplot(mesh,1);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
% end
% plot(pn(:,1),pn(:,2),'or','MarkerSize',7);
% plot(xu,yu,'-k','LineWidth',3);
% plot(xl,yl,'-k','LineWidth',3);
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_shockcurves.png
% 
% figure(2); clf; 
% simpplot(mesh.p, mesh.t);
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% print -dpng naca_mesh1.png
% 
% figure(3); clf; 
% simpplot(mesh2.p, mesh2.t);
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% print -dpng naca_mesh2.png

% figure(4); clf;
% scaplot(mesh,UDG1{3}(:,1,:),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% % xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% % ylabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_sol1.png
% 
% figure(5); clf;
% scaplot(mesh,ACG1{3}(:,1,:),[0 5e-3],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% % xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% % ylabel("$y$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; 
% colorbar off;
% colorbar('XTickLabel',{'0','0.001','0.002','0.003','0.004','0.005'}, ...
%    'XTick', 0:0.001:0.005, 'Location', 'northoutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_av1.png

% figure(6); clf;
% scaplot(mesh2,UDG2{6}(:,1,:),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_sol2.png

% figure(7); clf;
% scaplot(mesh2,ACG2{6}(:,1,:),[0 1.2e-4],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis tight; 
% colorbar off;
% colorbar('XTickLabel',{'0','0.00003','0.00006','0.00009','0.00012','0.00015'}, ...
%    'XTick', 0:0.00003:0.00015, 'Location', 'northoutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_av2.png

% figure(5); clf;
% scaplot(mesh,eulereval(UDG1{3},'p',1.4,0.85),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_pres1.png
% 
% figure(6); clf;
% scaplot(mesh2,eulereval(UDG2{6},'p',1.4,0.85),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% axis([-0.2 1.2 -0.2 0.8]);
% print -dpng naca_pres2.png

% figure(6); clf; hold on;
% plot(x, mach1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, mach2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, mach1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% plot(x, mach2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% text(0.65,1.2,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
% text(0.125,0.85,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$M$", 'interpreter', 'latex', 'FontSize', 22);
% box on;
% leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [25,10];
% print -dpng naca_mach3.png
% 
% figure(7); clf; hold on;
% plot(x, pres1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, pres2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, pres1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% plot(x, pres2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% text(0.67,0.8,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
% text(0.08,1.15,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$P$", 'interpreter', 'latex', 'FontSize', 22);
% box on;
% leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
% leg.ItemTokenSize = [25,10];
% print -dpng naca_pres3.png
% 
% figure(8); clf; hold on;
% plot(x, rho1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, rho2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
% plot(x, rho1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% plot(x, rho2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% text(0.67,0.8,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
% text(0.08,1.0,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 22);
% box on;
% leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
% leg.ItemTokenSize = [18,10];
% print -dpng naca_sol3.png
% 
% figure(1);clf; hold on;
% plot(xts(:,1),xts(:,2) ,'-', 'LineWidth', 2); 
% plot(x1,y1 ,'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
% leg = legend({'Exact shock curve', 'Predicted shock curve'}, 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [50,10];
% axis([0 0.7 0 1]);
% box on;
% print -dpng shockcurve.png
% 
% figure(2); clf;
% scaplot(mesh,UDG1(:,1,:),[0 2],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% print -dpng sol1.png
% 
% figure(3); clf;
% scaplot(mesh,ACG1(:,1,:),[0 0.026],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% print -dpng av1.png
% 
% figure(4); clf;
% scaplot(mesh2,UDG2(:,1,:),[0 2],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% print -dpng sol2.png
% 
% figure(5); clf;
% scaplot(mesh2,ACG2(:,1,:),[],2,0); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
% axis equal; axis tight; colorbar off; colorbar('NorthOutside');
% print -dpng av2.png
% 
% ii = [1 3 5 10];
% figure(6); clf; hold on;
% for i = 1:length(ii)
%     plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
%     plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
%     plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% end
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% text(0.065,1.825,'$\leftarrow t = 0.05$', 'interpreter', 'latex', 'FontSize',18);
% text(0.185,1.505,'$\leftarrow t = 0.20$', 'interpreter', 'latex', 'FontSize',18);
% text(0.325,1.275,'$\leftarrow t = 0.40$', 'interpreter', 'latex', 'FontSize',18);
% text(0.605,0.975,'$\leftarrow t = 0.90$', 'interpreter', 'latex', 'FontSize',18);
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
% ylabel("$u$", 'interpreter', 'latex', 'FontSize', 22);
% box on;
% leg = legend({'Exact solution', 'Initial numerical solution', 'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [40,10];
% print -dpng sol3.png
% 
% figure(7); clf; hold on;
% for i = 1:length(ii)
%     plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
%     plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
%     plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% end
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% box on;
% axis([-0.08 0.08 1.5 1.85]);
% print -dpng zoom1.png
% 
% figure(8); clf; hold on;
% for i = 1:length(ii)
%     plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
%     plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
%     plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% end
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% box on;
% axis([0.06 0.2 1.25 1.55]);
% print -dpng zoom2.png
% 
% figure(9); clf; hold on;
% for i = 1:length(ii)
%     plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
%     plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
%     plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% end
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% box on;
% axis([0.2 0.34 1.00 1.3]);
% print -dpng zoom3.png
% 
% figure(10); clf; hold on;
% for i = 1:length(ii)
%     plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
%     plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
%     plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% end
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% box on;
% axis([0.45 0.62 0.7 1.0]);
% print -dpng zoom4.png
% 
% figure(11); clf; 
% simpplot(mesh.p, mesh.t);
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% print -dpng mesh1.png
% 
% figure(12); clf; 
% simpplot(mesh2.p, mesh2.t);
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% print -dpng mesh2.png
