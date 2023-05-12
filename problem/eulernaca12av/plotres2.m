%load result1.mat

clear Cp1 x1;
for i=1:length(UDG1)
    mesh.dgnodes(:,3,:)=ACG1{i};
    [Cp1(:,i),Cf,x1,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh,app,UDG1{i},UH1{i},1);
end

clear Cp2 x2;
for i=1:length(UDG2)
    mesh2.dgnodes(:,3,:)=ACG2{i};
    [Cp2(:,i),Cf,x2,Cp2d,Cf2d,x2d,Ch,Ch2d]=getsurfacedata(master,mesh2,app,UDG2{i},UH2{i},1);
end

figure(1); clf; hold on;
plot(x1(:,1)/max(x1(:,1)), Cp1(:,4), '-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
plot(x2(:,1)/max(x2(:,1)), Cp2(:,8), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel("$x/c$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$C_p$", 'interpreter', 'latex', 'FontSize', 22);
axis tight; box on;
leg = legend({'Regular mesh',  'Shock-aligned mesh'}, 'FontSize', 16, 'Location', 'SE');
leg.ItemTokenSize = [25,10];
axis([0 1.02 -1.2 1.2]);
ax = gca;
fn = "naca_cp" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

ii = [1 3 5 7];
for i=1:length(ii)  
  figure(1); clf;
  scaplot(mesh,ACG1{ii(i)}(:,1,:),[],2,0); 
  set(gca,'FontSize',16);   
  axis equal; axis tight; axis off;
  colorbar off; colorbar('NorthOutside');
  axis([-0.1 1.1 -0.2 0.8]);
  ax = gca;
  fn = "naca_coarseav" + num2str(ii(i)) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end

ii = [1 3 5 7];
for i=1:length(ii)  
  figure(1); clf;
  scaplot(mesh,eulereval(UDG1{ii(i)},'p',1.4,0.85),[],2,0); 
  set(gca,'FontSize',16);   
  axis equal; axis tight; axis off;
  colorbar off; colorbar('NorthOutside');
  axis([-0.1 1.1 -0.2 0.8]);
  ax = gca;
  fn = "naca_coarsepressure" + num2str(ii(i)) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end

ii = [1 5 9 14];
ii=7;
for i=1:length(ii)  
  figure(1); clf;
  scaplot(mesh2,ACG2{ii(i)}(:,1,:),[],2,0); 
  set(gca,'FontSize',16);   
  axis equal; axis tight; axis off;
  colorbar off; colorbar('NorthOutside');
  axis([-0.1 1.1 -0.2 0.8]);
  ax = gca;
  fn = "naca_fineav" + num2str(ii(i)) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end

ii=14;
for i=1:length(ii)  
  figure(1); clf;
  scaplot(mesh2,eulereval(UDG2{ii(i)},'M',1.4,0.85),[],2,0); 
  set(gca,'FontSize',16);   
  axis equal; axis tight; axis off;
  colorbar off; colorbar('NorthOutside');
  axis([-0.1 1.1 -0.2 0.8]);
  ax = gca;
  fn = "naca_finemach" + num2str(ii(i)) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end


% figure(2); clf; 
% meshplot(mesh,1);
% set(gca,'FontSize',18); 
% set(gca,'LooseInset',get(gca,'TightInset'));
% axis on; box on; axis tight;
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
% ylabel("$y$", 'interpreter', 'latex', 'FontSize', 24);
% ax = gca;
% fn = "naca_mesh" + num2str(1) + ".png";
% exportgraphics(ax,fn,'Resolution',200); 
% 
% figure(2); clf; 
% meshplot(mesh2,1);
% set(gca,'FontSize',18); 
% set(gca,'LooseInset',get(gca,'TightInset'));
% xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
% ylabel("$y$", 'interpreter', 'latex', 'FontSize', 24);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% ax = gca;
% fn = "naca_mesh" + num2str(2) + ".png";
% exportgraphics(ax,fn,'Resolution',200); 

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
