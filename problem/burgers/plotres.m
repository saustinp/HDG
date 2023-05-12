load results.mat

figure(1);clf; hold on;
plot(xts(:,1),xts(:,2) ,'-', 'LineWidth', 2); 
plot(x1,y1 ,'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
leg = legend({'Exact shock curve', 'Predicted shock curve'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [50,10];
axis([0 0.7 0 1]);
box on;
print -dpng shockcurve.png

figure(2); clf;
scaplot(mesh,UDG1(:,1,:),[0 2],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
axis equal; axis tight; colorbar off; colorbar('NorthOutside');
print -dpng sol1.png

figure(3); clf;
scaplot(mesh,ACG1(:,1,:),[0 0.026],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
axis equal; axis tight; colorbar off; colorbar('NorthOutside');
print -dpng av1.png

figure(4); clf;
scaplot(mesh2,UDG2(:,1,:),[0 2],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
axis equal; axis tight; colorbar off; colorbar('NorthOutside');
print -dpng sol2.png

figure(5); clf;
scaplot(mesh2,ACG2(:,1,:),[],2,0); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
axis equal; axis tight; colorbar off; colorbar('NorthOutside');
print -dpng av2.png

ii = [1 3 5 10];
figure(6); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
text(0.065,1.825,'$\leftarrow t = 0.05$', 'interpreter', 'latex', 'FontSize',18);
text(0.185,1.505,'$\leftarrow t = 0.20$', 'interpreter', 'latex', 'FontSize',18);
text(0.325,1.275,'$\leftarrow t = 0.40$', 'interpreter', 'latex', 'FontSize',18);
text(0.605,0.975,'$\leftarrow t = 0.90$', 'interpreter', 'latex', 'FontSize',18);
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$u$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Exact solution', 'Initial numerical solution', 'Final numerical solution'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [40,10];
print -dpng sol3.png

figure(7); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
box on;
axis([-0.08 0.08 1.5 1.85]);
print -dpng zoom1.png

figure(8); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
box on;
axis([0.06 0.2 1.25 1.55]);
print -dpng zoom2.png

figure(9); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
box on;
axis([0.2 0.34 1.00 1.3]);
print -dpng zoom3.png

figure(10); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux1(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux2(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
box on;
axis([0.45 0.62 0.7 1.0]);
print -dpng zoom4.png

figure(11); clf; 
simpplot(mesh.p, mesh.t);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
print -dpng mesh1.png

figure(12); clf; 
simpplot(mesh2.p, mesh2.t);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
print -dpng mesh2.png

figure(12); clf; 
meshplot(mesh2, 1);
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'));
axis on; axis equal; axis tight;
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
ax = gca;
fn = "bg_mesh2" + ".png";
exportgraphics(ax,fn,'Resolution',200); 

for i = 1:length(ACG2)
  figure(5); clf;
  scaplot(mesh2,ACG2{i}(:,1,:),[],2,0); 
  set(gca,'FontSize',16); 
  xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
  ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
  axis equal; axis tight; colorbar off; colorbar('NorthOutside');
  ax = gca;
  fn = "bg_av" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end

for i = 1:length(ACG2)
  figure(5); clf;
  scaplot(mesh2,UDG2{i}(:,1,:),[0 2],2,0); 
  set(gca,'FontSize',16); 
  xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
  ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
  axis equal; axis tight; colorbar off; colorbar('NorthOutside');
  ax = gca;
  fn = "bg_sol" + num2str(i) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 
end


  figure(5); clf;
  scaplot(mesh,ACG1{4}(:,1,:),[],2,0); 
  set(gca,'FontSize',16); 
  xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
  ylabel("$t$", 'interpreter', 'latex', 'FontSize', 22);
  axis equal; axis tight; colorbar off; colorbar('NorthOutside');
  ax = gca;
  fn = "bg_coarsesav" + num2str(4) + ".png";
  exportgraphics(ax,fn,'Resolution',200); 


