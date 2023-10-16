    %% NOTE: NOT COMPLETED YET
    ncu = app.ncu;
    ne = mesh.ne;
    npe = master.npe;
    
    [pg, Xx, jac] = volgeom(master.shapmv,permute(mesh.dgnodes,[1 3 2]));
    M = app.fc_q*reshape(master.shapvgdotshapvl(:,:,1)*reshape(jac,[master.ngv mesh.ne]),[master.npv master.npv mesh.ne]);
    
    i_train = 2:2:length(Ma_array);
    Ma_array_train = Ma_array(i_train);
    p_rb = length(Ma_array_train)

    phi_Ma={};
    for i = 1:p_rb
        str = "mach_sweep_data/Mach"+string(Ma_array_train(i));
        wkspc_tmp = load(str+"_res");
        phi_Ma{i} = wkspc_tmp.q1;
    end

    Phi_x = zeros(npe*ne,p_rb); 
    MPhi_x = Phi_x;
    Phi_y = zeros(npe*ne,p_rb);
    MPhi_y = Phi_y;
    Phi_UDG1 = zeros(npe*1*ne,p_rb); 
    Phi_UDG2_smooth = zeros(npe*1*ne,p_rb);
    Phi_UDG2_sharp = zeros(npe*1*ne,p_rb);
    Phi_UDG3_smooth = zeros(npe*1*ne,p_rb);
    Phi_UDG3_sharp = zeros(npe*1*ne,p_rb);
    for i = 1:p_rb
        q = phi_Ma{i};
        UDG1tmp = UDG1struct{i_train(i)};
        UDG1tmp = UDG1tmp{end}(:,1,:);
        UDG2tmp = UDG2struct{i_train(i)};
%         UDG3tmp = UDG3struct_fix{i_train(i)};

        UDG2tmp_smooth = UDG2tmp{1}(:,1,:);
        UDG2tmp_sharp = UDG2tmp{4}(:,1,:);

%         UDG2tmp_sharp = UDG2tmp(:,1,:);

        Phi_UDG1(:,i) = UDG1tmp(:);
        
        Phi_UDG2_smooth(:,i) = UDG2tmp_smooth(:);
%         MPhi_UDG2_smooth(:,i) = apply_weight_to_vector(M, Phi_UDG2_smooth(:,i),mesh,master);
        Phi_UDG2_sharp(:,i) = UDG2tmp_sharp(:);
%         MPhi_UDG2_sharp(:,i) = apply_weight_to_vector(M, Phi_UDG2_sharp(:,i),mesh,master);
       

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
xlim([1 p_rb])
title("$\Phi_x$","Interpreter","latex"); 
set(gca,'FontSize',15); grid on
subplot(1,2,2); semilogy(svals_y/svals_y(1), 'LineWidth',2); 
xlim([1 p_rb])
title("$\Phi_y$","Interpreter","latex");
set(gca,'FontSize',15); grid on

%%
[~,svals_mapped_smooth2] = snapshot_POD(Phi_UDG2_smooth);
[Wmodes2,svals_mapped_sharp2] = snapshot_POD(Phi_UDG2_sharp);
% 
% [Wmodes3_smooth,svals_mapped_smooth3] = snapshot_POD(Phi_UDG3_smooth);
% [Wmodes3_sharp,svals_mapped_sharp3] = snapshot_POD(Phi_UDG3_sharp);

[Wmodes1,svals_fixed] = snapshot_POD(Phi_UDG1);
figure(3); semilogy(svals_fixed/svals_fixed(1), 'LineWidth',2)
figure(3); hold on; semilogy(svals_mapped_sharp2/svals_mapped_sharp2(1), 'LineWidth',2)
% figure(3); hold on; semilogy(svals_mapped_smooth2/svals_mapped_smooth2(1), 'LineWidth',2)
% figure(3); hold on; semilogy(svals_mapped_sharp3/svals_mapped_sharp3(1), 'LineWidth',2)
% figure(3); hold on; semilogy(svals_mapped_smooth3/svals_mapped_smooth3(1), 'LineWidth',2)
title("Singular value decay: cylinder Ma 3 to 11")
legend(["Fixed mesh", "Adaptive Mesh 1", "Adaptive Mesh 2"])
set(gca,'FontSize',15); grid on

%%
%%
for i = 1:p_rb
    tst1 = Wmodes1(:,i);
%     tst1 = reshape(tst1, [npe ncu ne]);
    tst1 = reshape(tst1, [npe 1 ne]);

    tst2 = Wmodes3_sharp(:,i);
%     tst2 = reshape(tst2, [npe ncu ne]);
        tst2 = reshape(tst2, [npe 1 ne]);

    figure(4); clf;
    subplot(1,2,1);
    scaplot(mesh, tst1,[],2);

% scaplot(mesh, eulereval(tst1,'p',1.4,Minf),[-2 2]);
    subplot(1,2,2); 
%     scaplot(mesh, eulereval(tst2,'p',1.4,Minf),[],2);
    scaplot(mesh, tst2,[],2);
    drawnow; waitforbuttonpress
end



%%
Ma_array_test = Ma_array;
% Ma_array_test = [3.5, 5, 7, 11]
for i = 1:length(Ma_array_test)
    str = "mach_sweep_data/Mach"+string(Ma_array_test(i));
    disp("Ma = " + string(Ma_array_test(i)))
    wkspc_tmp = load(str+"_res");
    
    UDG1_tmp = wkspc_tmp.UDG1{end}(:,1,:);
%     figure(1); scaplot(mesh, UDG1_tmp);
%     tst = Wmodes1*(pinv(Wmodes1)*UDG1tmp(:));
    alpha = pinv(Wmodes1)*UDG1_tmp(:);
%     disp(alpha);
    tst = Wmodes1*(alpha);
    figure(1); clf; scaplot(mesh, tst)
        figure(2); clf; scaplot(mesh, abs(tst(:) - UDG1_tmp(:)))
%     disp(norm(tst(:)-UDG1_tmp(:)))
    [err, err_rel] = calerror_ROM(reshape(tst, size(UDG1_tmp)), UDG1_tmp, mesh, master);
    disp("(E, E_r) = " + string(err) +", " + string(err_rel) + ")")
    drawnow; waitforbuttonpress;
end

%%
% √=Ma_array_test = Ma_array;
for i = 1:length(Ma_array_test)
    str = "mach_sweep_data/Mach"+string(Ma_array_test(i));
    disp("Ma = " + string(Ma_array_test(i)))
    wkspc_tmp = load(str+"_res");
%     mesh_tmp = mesh;
    UDG1_tmp = wkspc_tmp.UDG2{4}(:,1,:);
    mesh_tmp = wkspc_tmp.mesh2;
%     figure(1); scaplot(mesh, UDG1_tmp);
%     tst = Wmodes1*(pinv(Wmodes1)*UDG1tmp(:));
    alpha = pinv(Wmodes2)*UDG1_tmp(:);
%     disp(alpha);
    tst = Wmodes2*(alpha);
    figure(1); clf; scaplot(mesh_tmp, tst)
    figure(2); clf; scaplot(mesh_tmp, abs(tst(:) - UDG1_tmp(:)))
    [err, err_rel] = calerror_ROM(UDG1_tmp, reshape(tst, size(UDG1_tmp)), mesh_tmp, master);
    disp("(E, E_r) = " + string(err) +", " + string(err_rel) + ")")
% 6 2.5 1 1.5
% 1.5 0.5 0.5 0.5
    drawnow; waitforbuttonpress;
end