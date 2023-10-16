porder = 4;
% 
% meshstruct = {};
% UDG1struct = {};
% UDG2struct = {};
iter = 1;
Ma_array = [2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19];
for Mach = Ma_array(iter:end)
    disp("========Running Mach number=====" + string(Mach));
    str = "mach_sweep_data/Mach"+string(Mach);
    nstage = 1;
    torder = 1;
    % Mach   = 5;
    aoa    = 0.0;
    hybrid = 'hdg';
    
    gam = 1.4;
    epslm = 0.0;
    Minf = Mach;                  % Infinity conditions
    pinf = 1/(gam*Minf^2);
    Re = inf;
    Pr = 0.72;
    alpha = aoa*pi/180;
    tau = 1;
    
    nd    = 2;
    ntime = 30;
    dt = 1e-3*2.^(0:ntime);
    dt = repmat(dt,[nd 1]);
    dt = dt(:);
    
    ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
    
    clear app;
    app.source = 'source';
    app.flux = 'flux';
    app.fbou = 'fbou';
    app.fhat = 'fhat';
    app.hybrid = hybrid;
    app.localsolve = 1;
    app.arg = {gam,Minf,epslm,tau};
    app.bcm  = [5,2,5,6];  
    app.bcs  = [ui; ui; ui; ui];
    app.bcd  = [1,1,1,1];  
    app.bcv  = [0; 0; 0; 0];
    
    app.tdep = false;
    app.wave = false;
    app.alag = false;
    app.flg_q = 1;
    app.flg_p = 0;
    app.flg_g = 0;
    
    app.ndim = 2;
    app.nch  = 2+app.ndim;                % Number of componets of UH
    app.nc   = app.nch*3;                   % Number of componeents of UDG
    app.ncu  = app.nch;                   % Number of components of U
    
    app.time = [];
    app.dtfc = [];
    app.alpha = [];
    
    app.fc_q = 1;
    app.fc_u = 0;
    app.tdep = false;
    app.adjoint = 0;
    
    meshsq = mkmesh_square(41,41,porder,1,1,1,1,1);
    meshsq.p(:,1) = logdec(meshsq.p(:,1),0.5);
    meshsq.dgnodes(:,1,:) = logdec(meshsq.dgnodes(:,1,:),0.5);
    mesh = mkmesh_halfcircle(meshsq,1,4,7,pi/2,3*pi/2);
    
    master = mkmaster(mesh,2*porder);
    [master,mesh] = preprocess(master,mesh,hybrid);
    mesh = mkcgmesh(mesh);
    mesh.ib = [];
    mesh.in = 1:size(mesh.p2,1);
    x = mesh.p2(:,1); y = mesh.p2(:,2);
    mesh.ib = find(abs(x.^2 + y.^2)<1+1e-5);
    mesh.in = setdiff( (1:size(mesh.p2,1))', mesh.ib);
    [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
    mesh.dist = tanh(meshdist(mesh,2)*25);
    
    UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
    UDG(:,2,:) = UDG(:,2,:).*tanh(meshdist(mesh,2)*4);
    UDG(:,3,:) = UDG(:,3,:).*tanh(meshdist(mesh,2)*4);
    TnearWall = pinf/(gam-1); % Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
    UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
    UH = inituhat(master,mesh.elcon,UDG,app.ncu);
    SH = [];
    
    mesh.dgnodes(:,3,:) = 0.06*mesh.dist;
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
    
    mesh.dgnodes(:,3,:) = 0.04*mesh.dist;
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
    
    mesh.dgnodes(:,3,:) = 0.03*mesh.dist;
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
    
    mesh.dgnodes(:,3,:) = 0.03*tanh(meshdist(mesh,2)*5);
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
    
    app.S0=0.2; app.lambda = 0.03; app.kappa=3;
    a = avf(mesh, master, app, UDG);
    figure(1); clf; scaplot(mesh, a,[],2); 
    % Initial solve
    [UDG1, UH1, ACG1, mine1, minf1, ming1, mesh1] = hdgsolve_avloop2(master, mesh, app, mesh.dgnodes(:,1:2,:), UDG, UH, 0.04, 4, 3);
    figure(1); clf; scaplot(mesh, eulereval(UDG1{2},'M',gam,Minf),[],2); 
    
    % Solve Monge Ampere
    [u1,q1,uhat1,v1,rho1,iter1] = hdgsolve_mah(master, mesh, mesh1, UDG1{end}, 20e-3, []);
    [elist, xi] = locatexinmesh(mesh1, q1, [], 1e-4);
    % Map to new mesh
    UDG0 = evalfield(mesh1, UDG1{end}, elist, xi);
    UDG0 = permute(reshape(UDG0, [master.npe mesh.ne 12]),[1 3 2]); 
    UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);
    % Solve on new mesh
    [UDG2, UH2, ACG2, mine2, minf2, ming2, mesh2] = hdgsolve_avloop(master, mesh, app, q1, UDG0, UH0, 0.02, 5);

    meshstruct{iter} = mesh2;
    UDG1struct{iter} = UDG1;
    UDG2struct{iter} = UDG2;
    iter = iter + 1;

    save(str+"_res");

%     figure(3); clf; scaplot(mesh2, eulereval(UDG2{end},'M',1.4,Minf),[],1,2); axis on; axis equal; axis tight;
%     mf1 = 18;
%     set(gca,'FontSize',mf1); 
%     box on; axis tight;
%     ax = gca;
%     fn = str + "_mapped_mesh" + ".png";
%     exportgraphics(ax,fn,'Resolution',200); 
%     figure(3); clf; scaplot(mesh, eulereval(UDG2{2},'M',1.4,Minf),[]); axis on; axis equal; axis tight;
%     set(gca,'FontSize',mf1); 
%     box on; axis tight;
%     ax = gca;
%     fn = str + "_fixed_mesh" + ".png";
%     exportgraphics(ax,fn,'Resolution',200); 
%     waitforbuttonpress();
end
    % mesht = mesh; mesht.dgnodes = q; figure(1); clf; meshplot(mesht,1); axis on; axis equal; axis tight;
    % mesht = mesh; mesht.dgnodes = q1; figure(2); clf; meshplot(mesht,1); axis on; axis equal; axis tight;
    % figure(3); clf; scaplot(mesh1, eulereval(UDG1{5},'M',1.4,7),[],1);  axis on; axis equal; axis tight;
    % figure(4); clf; scaplot(mesht, eulereval(UDG0,'M',1.4,7),[],1);  axis on; axis equal; axis tight;
    % [u2,q2,uhat2,v2,rho2,iter2] = hdgsolve_mah(master, mesh, mesh1, UDG1{4}, 30e-3, [], 0.25);
    % [elist, xi] = locatexinmesh(mesh1, q2, [], 1e-4);
    % UDG0 = evalfield(mesh1, UDG1{4}, elist, xi);
    % UDG0 = permute(reshape(UDG0, [master.npe mesh.ne 12]),[1 3 2]);
    % UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);
    % [UDG3, UH3, ACG3, mine3, minf3, ming3, mesh3] = hdgsolve_avloop(master, mesh, app, q2, UDG0, UH0, 0.02, 5);
    % % 

    %% Quick plotting check
%     Ma_array = [11, 15];
    pinf_array = 1.0 ./ (1.4*Ma_array.^2)
    Tinf_array = pinf_array / (gam-1);
%     UDGma15 = UDG2struct{end};
%     UDGma7 = UDG2struct{1};
%     UDGtmpstruct = {UDGma7{end}, UDGma15{end}};
%     meshtmpstruct = {meshstruct{1}, meshstruct{end}};

    x1 = cell(2,1);
    r1 = cell(2,1);
    p1 = cell(2,1);
    m1 = cell(2,1);
    t1 = cell(2,1);
    for i = 1:7
        Minf = Ma_array(i);
        UDGtmp_struct = UDG2struct{i};
        meshtmp = meshstruct{i};
        UDGtmp = UDGtmp_struct{end};
%         meshtmp = mesh;
      [~,r1{i}] = getfieldaty(meshtmp,eulereval(UDGtmp,'r',1.4,Minf));    
      [~,p1{i}] = getfieldaty(meshtmp,eulereval(UDGtmp,'p',1.4,Minf));    
      [x1{i},m1{i}] = getfieldaty(meshtmp,eulereval(UDGtmp,'M',1.4,Minf));    
      [~,t1{i}] = getfieldaty(meshtmp,eulereval(UDGtmp,'T',1.4,Minf));   

      [~,r1_refmesh{i}] = getfieldaty(mesh,eulereval(UDGtmp,'r',1.4,Minf));    
      [~,p1_refmesh{i}] = getfieldaty(mesh,eulereval(UDGtmp,'p',1.4,Minf));    
      [x1_refmesh{i},m1_refmesh{i}] = getfieldaty(mesh,eulereval(UDGtmp,'M',1.4,Minf));    
      [~,t1_refmesh{i}] = getfieldaty(mesh,eulereval(UDGtmp,'T',1.4,Minf));  
    end

    
    % gam = 1.4;
    % Mach = 17.605;
    % Minf = Mach;
    % pinf  = 1/(gam*Minf^2);
    % Tinf = pinf/(gam-1);
    % rinf = 1.0;
    
    figure(20); clf; hold on;
    for i = 1:7
%          UDG2struct{1}
            plot(x1{i},r1{i},'LineWidth',2,'MarkerSize',8);
            xlabel("$x$", 'Interpreter','latex')
            ylabel("\rho")
            set(gca,'FontSize',18); 
            leg = legend(["$M_{\infty} = 2$", "$M_{\infty} = 3$", "$M_{\infty} = 4$", "$M_{\infty} = 5$", "$M_{\infty} = 7$", "$M_{\infty} = 9$", "$M_{\infty} = 11$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
            leg.ItemTokenSize = [20,10];
    end

    figure(30); clf; hold on;
    for i = 1:7
%          UDG2struct{1}
            plot(x1_refmesh{i},r1_refmesh{i},'LineWidth',2,'MarkerSize',8);
            xlabel("$\mathcal{G}(x)$", 'Interpreter','latex')
            ylabel("\rho")
            set(gca,'FontSize',18); 
            leg = legend(["$M_{\infty} = 2$", "$M_{\infty} = 3$", "$M_{\infty} = 4$", "$M_{\infty} = 5$", "$M_{\infty} = 7$", "$M_{\infty} = 9$", "$M_{\infty} = 11$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
            leg.ItemTokenSize = [20,10];
%             axis([-3 -1 1 160]);
    end

    %%
    cnt = 1;
    for i = [2 4 7]
%          UDG2struct{1}
%         figure(i); clf; 
        Minf = Ma_array(i);
        UDGtmp_struct = UDG2struct{i};
        meshtmp = meshstruct{i};
        UDGtmp = UDGtmp_struct{1};
        figure(1); subplot(1,3,cnt); scaplot(mesh, eulereval(UDGtmp,'M',1.4,Minf), [0 Minf]); colormap('jet')
        title('$\mu=$'+string(Minf),'Interpreter','latex')
        set(gca,'FontSize',18); 
        axis on; box on; axis tight;

        figure(2); subplot(1,3,cnt); scaplot(meshtmp, eulereval(UDGtmp,'M',1.4,Minf),[0 Minf]); colormap('jet')
%         title("Mach " + string(Ma_array(i)))
%         set(gca,'FontSize',18);
        title('$\mu=$'+string(Minf),'Interpreter','latex')
        set(gca,'FontSize',18); 
        axis on; box on; axis tight;

        figure(3); subplot(1,3,cnt); meshplot(meshtmp,1)
%         title("Mach " + string(Ma_array(i)))
%         set(gca,'FontSize',18);
        title('$\mu=$'+string(Minf),'Interpreter','latex')
        set(gca,'FontSize',18); 
        axis on; box on; axis tight;


        cnt = cnt + 1;
%             plot(x1_refmesh{i},r1_refmesh{i},'LineWidth',2,'MarkerSize',8);
%             xlabel("$\mathcal{G}(x)$", 'Interpreter','latex')
%             ylabel("\rho")
%             set(gca,'FontSize',18); 
%             leg = legend(["$M_{\infty} = 2$", "$M_{\infty} = 3$", "$M_{\infty} = 4$", "$M_{\infty} = 5$", "$M_{\infty} = 7$", "$M_{\infty} = 9$", "$M_{\infty} = 11$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
%             leg.ItemTokenSize = [20,10];
%             axis([-3 -1 1 160]);
    end


%     plot(x1{1},r1{1},':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
%     plot(x1{2},r1{2},'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    % plot(x1{5},r1{5},'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
    % plot(x1{7},r1{7},'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
    % set(gca,'FontSize',18); 
    % axis tight; axis on; box on; 
    % xlabel("$x$", 'interpreter', 'latex', 'FontSize', 20);
    % ylabel("$\rho/\rho_\infty$", 'interpreter', 'latex', 'FontSize', 20);
    % leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
    % leg.ItemTokenSize = [30,10];
    % axis([-3 -1 1 160]);
    % ax = gca;
    % exportgraphics(ax,"nscyl18_adaptive1_centerlinedensity.png",'Resolution',200); 
    
%     figure(30); clf; hold on;
%     plot(x1{1},p1{1}/pinf_array(1),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
%     plot(x1{2},p1{2}/pinf_array(2),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    
    
%     figure(40); clf; hold on;
%     plot(x1{1},t1{1}/Tinf_array(1),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
%     plot(x1{2},t1{2}/Tinf_array(2),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    % plot(x1{5},t1{5}/Tinf,'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
    % plot(x1{7},t1{7}/Tinf,'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
    % set(gca,'FontSize',18); 
    % axis tight; axis on; box on; 
    % xlabel("$x$", 'interpreter', 'latex', 'FontSize', 20);
    % ylabel("$T/T_\infty$", 'interpreter', 'latex', 'FontSize', 20);
    % leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
    % leg.ItemTokenSize = [30,10];
    % axis([-3 -1 [0.005761 0.38]/Tinf]);
    % ax = gca;
    % exportgraphics(ax,"nscyl18_adaptive1_centerlinetemperature.png",'Resolution',200); 
    
    % ExpSt = csvread('compar_data/LAURA_St.txt');
    % ExpCf = csvread('compar_data/LAURA_Cf.txt');
    % ExpCp = csvread('compar_data/LAURA_Cp.txt');
    
    % wid=2;
    % gam1=gam-1;
    % Tinf = pinf/(gam-1);
    % Ttinf = Tinf * (1 + (gam-1)/2 * Minf^2);
    % TisoW = 500./200. * Tinf;
    % deltaT = Ttinf - TisoW;
    % elemAvg = 0;
    % Cp1 = cell(7,1);
    % Cf1 = cell(7,1);
    % Ch1 = cell(7,1);
    % x1 = cell(7,1);
    % theta1 = cell(7,1);
    % for i = 1:7
    %   [Cp1{i},Cf1{i},x1{i},~,~,~,Ch1{i}]=getsurfacedata(master,mesh,app,UDG1{i},UH1{i},wid,elemAvg,deltaT);
    %   theta1{i} = atan2(x1{i}(:,2),-x1{i}(:,1));
    %   ii = theta1{i}>0; Cf1{i}(ii) = -Cf1{i}(ii);  
    % end
    
    % figure(2); clf; hold on;
    % plot(theta1{1}*180/pi,-Cp1{1},':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{3}*180/pi,-Cp1{3},'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{5}*180/pi,-Cp1{5},'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
    % plot(theta1{7}*180/pi,-Cp1{7},'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
    % plot(ExpCp(:,1),ExpCp(:,2),'k-','LineWidth',2,'MarkerSize',8);
    % set(gca,'FontSize',18); 
    % axis tight; axis on; box on; 
    % xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
    % ylabel("$C_p$", 'interpreter', 'latex', 'FontSize', 20);
    % leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{LAURA}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
    % leg.ItemTokenSize = [30,10];
    % axis([-90 90 0 1.9]);
    % ax = gca;
    % xticks([-90:30:90]);
    % yticks([0:0.3:1.8]);
    % exportgraphics(ax,"nscyl18_adaptive1_pressurecoefficient.png",'Resolution',200); 
    
    % figure(2); clf; hold on;
    % plot(theta1{1}*180/pi,Ch1{1},':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{3}*180/pi,Ch1{3},'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{5}*180/pi,Ch1{5},'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
    % plot(theta1{7}*180/pi,Ch1{7},'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
    % plot(ExpSt(:,1),ExpSt(:,2),'k-','LineWidth',2,'MarkerSize',8);
    % set(gca,'FontSize',18); 
    % axis tight; axis on; box on; 
    % xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
    % ylabel("$C_h$", 'interpreter', 'latex', 'FontSize', 20);
    % leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{LAURA}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
    % leg.ItemTokenSize = [30,10];
    % axis([-90 90 0 9e-3]);
    % ax = gca;
    % xticks([-90:30:90]);
    % yticks([0:1e-3:9e-3]);
    % exportgraphics(ax,"nscyl18_adaptive1_stanton.png",'Resolution',200); 
    
    % figure(2); clf; hold on;
    % plot(theta1{1}*180/pi,Cf1{1},':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{3}*180/pi,Cf1{3},'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
    % plot(theta1{5}*180/pi,Cf1{5},'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
    % plot(theta1{7}*180/pi,Cf1{7},'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
    % plot(ExpCf(:,1),ExpCf(:,2),'k-','LineWidth',2,'MarkerSize',8);
    % set(gca,'FontSize',18); 
    % axis tight; axis on; box on; 
    % xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
    % ylabel("$C_f$", 'interpreter', 'latex', 'FontSize', 20);
    % leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{LAURA}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
    % leg.ItemTokenSize = [30,10];
    % axis([-90 90 -6e-3 6e-3]);
    % ax = gca;
    % xticks([-90:30:90]);
    % yticks([-6e-3:2e-3:6e-3]);
    % exportgraphics(ax,"nscyl18_adaptive1_skinfriction.png",'Resolution',200); 
    %% QUick POD check of MA outputs
    ncu = app.ncu;
    ne = mesh.ne;
    npe = master.npe;
    p_rb = 7;
    [pg, Xx, jac] = volgeom(master.shapmv,permute(mesh.dgnodes,[1 3 2]));
    M = app.fc_q*reshape(master.shapvgdotshapvl(:,:,1)*reshape(jac,[master.ngv mesh.ne]),[master.npv master.npv mesh.ne]);

    phi_Ma={};
    for i = 1:p_rb
        str = "mach_sweep_data/Mach"+string(Ma_array(i));
        wkspc_tmp = load(str+"_res");
        phi_Ma{i} = wkspc_tmp.q1;
    end

    Phi_x = zeros(npe*ne,p_rb); 
    MPhi_x = Phi_x;
    Phi_y = zeros(npe*ne,p_rb);
    MPhi_y = Phi_y;
    Phi_UDG1 = zeros(npe*ncu*ne,p_rb); 
    Phi_UDG2_smooth = zeros(npe*1*ne,p_rb);
    Phi_UDG2_sharp = zeros(npe*1*ne,p_rb);
    for i = 1:p_rb
        q = phi_Ma{i};
        UDG1tmp = UDG1struct{i};
        UDG1tmp = UDG1tmp{end}(:,1:ncu,:);
        UDG2tmp = UDG2struct{i};
        UDG2tmp_smooth = UDG2tmp{1}(:,4,:);
        UDG2tmp_sharp = UDG2tmp{end}(:,4,:);
%         UDG2tmp_sharp = UDG2tmp(:,1,:);

        Phi_UDG1(:,i) = UDG1tmp(:);
        Phi_UDG2_smooth(:,i) = UDG2tmp_smooth(:);
        MPhi_UDG2_smooth(:,i) = apply_weight_to_vector(M, Phi_UDG2_smooth(:,i),mesh,master);
        Phi_UDG2_sharp(:,i) = UDG2tmp_sharp(:);
        MPhi_UDG2_sharp(:,i) = apply_weight_to_vector(M, Phi_UDG2_sharp(:,i),mesh,master);


        qx = q(:,1,:);
        qy = q(:,2,:);
        qx = qx(:);
        qy = qy(:);
        Phi_x(:,i) = qx(:);
        MPhi_x(:,i) = apply_weight_to_vector(M, Phi_x(:,i),mesh,master);
        Phi_y(:,i) = qy(:);
        MPhi_y(:,i) = apply_weight_to_vector(M, Phi_y(:,i),mesh,master);
    end

[V_x,D_x] = eig(Phi_x'*MPhi_x); %%TODO: should multiply by mass matrix to have correct inner product
[~,ii] = sort(diag(D_x),'descend');
D_x = D_x(ii,ii);
V_x = V_x(:,ii);
Wmodes_x = Phi_x*V_x;
svals_x = diag(D_x);
    
[V_y,D_y] = eig(Phi_y'*MPhi_y);
[~,jj] = sort(diag(D_y),'descend');
D_y = D_y(jj,jj);
V_y = V_y(:,jj);
Wmodes_y = Phi_y*V_y;
svals_y = diag(D_y);
%%
figure(100); clf; 
subplot(1,2,1); semilogy(svals_x/svals_x(1), 'LineWidth',2); 
xlim([1 7])
title("$\Phi_x$","Interpreter","latex"); 
set(gca,'FontSize',15); grid on
subplot(1,2,2); semilogy(svals_y/svals_y(1), 'LineWidth',2); 
xlim([1 7])
title("$\Phi_y$","Interpreter","latex");
set(gca,'FontSize',15); grid on

%%
[Wmodes2,svals_mapped_smooth] = snapshot_POD(Phi_UDG2_smooth,M,mesh,master);
[~,svals_mapped_sharp] = snapshot_POD(Phi_UDG2_sharp,M,mesh,master);

[Wmodes1,svals_fixed] = snapshot_POD(Phi_UDG1,M,mesh,master);
figure(3); semilogy(svals_fixed/svals_fixed(1), 'LineWidth',2)
figure(3); hold on; semilogy(svals_mapped_sharp/svals_mapped_sharp(1), 'LineWidth',2)
figure(3); hold on; semilogy(svals_mapped_smooth/svals_mapped_smooth(1), 'LineWidth',2)

%%
for i = 1:7
    tst1 = Wmodes1(:,i);
    tst1 = reshape(tst1, [npe ncu ne]);
%     tst1 = reshape(tst1, [npe 1 ne]);

    tst2 = Wmodes2(:,i);
%     tst2 = reshape(tst2, [npe ncu ne]);
        tst2 = reshape(tst2, [npe 1 ne]);

    figure(4); clf;
    subplot(1,2,1); 
%     scaplot(mesh, tst1(:,4,:),[],2);

scaplot(mesh, eulereval(tst1,'p',1.4,Minf),[-2 2]);
    subplot(1,2,2); 
%     scaplot(mesh, eulereval(tst2,'p',1.4,Minf),[],2);
    scaplot(mesh, tst2,[],2);
    drawnow; waitforbuttonpress
end



