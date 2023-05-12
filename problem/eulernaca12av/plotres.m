%load results.mat

yu = linspace(0.04326, 0.7172, 101);
xu = polyval(poly,yu);

yl = linspace(-0.16794,-0.0597, 101);
xl = polyval(poly2,yl);

figure(1); clf; hold on;
meshplot(mesh,1);
for i = 1:size(pp,1)
    plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
end
plot(pn(:,1),pn(:,2),'or','MarkerSize',7);
plot(xu,yu,'-k','LineWidth',3);
plot(xl,yl,'-k','LineWidth',3);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight;
axis([-0.2 1.2 -0.2 0.8]);
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

figure(6); clf; hold on;
plot(x, mach1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, mach2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, mach1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
plot(x, mach2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
text(0.65,1.2,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
text(0.125,0.85,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\mbox{Mach number}$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Regular mesh',  'Shock-aligned mesh'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [25,10];
print -dpng naca_mach3.png

figure(7); clf; hold on;
plot(x, pres1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, pres2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, pres1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
plot(x, pres2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
text(0.67,0.8,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
text(0.08,1.15,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\mbox{Pressure}$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Regular mesh',  'Shock-aligned mesh'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [25,10];
print -dpng naca_pres3.png

figure(8); clf; hold on;
plot(x, rho1(:,8), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, rho2(:,8), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x, rho1(:,7), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
plot(x, rho2(:,7), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
text(0.67,0.8,'$\leftarrow y = 0.15$', 'interpreter', 'latex', 'FontSize',18);
text(0.08,1.0,'$\leftarrow y = -0.15$', 'interpreter', 'latex', 'FontSize',18);
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$\rho$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Initial numerical solution',  'Final numerical solution'}, 'FontSize', 16, 'Location', 'SW');
leg.ItemTokenSize = [18,10];
print -dpng naca_sol3.png

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
