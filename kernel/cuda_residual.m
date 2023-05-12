function Ru = cuda_residual(Mi, Ru, q, Rq, udg, sdg, xdg, uh, uinf, tempen, tempeg, tempfn, tempfg, shapen, shapeg, shapfn, shapfg, param, time, fc_q, fc_u, tdep, f2e, nn, nme, nmf)
% master: shapen, shapeg, shapfn, shapfg 
% mesh: Mi, xdg, f2e, nme, nmf 
% sol: udg, sdg, odg, uh, Ru, Rq
% temp: q, tempen, tempeg, tempfn, tempfg
% app: param, time, fc_q, fc_u, tdep, uinf, nn

ncu = nn(3);
q = cuda_getq(q, Mi, Rq, sdg, udg, xdg, uh, uinf, tempen, tempeg, tempfn, tempfg, shapen, shapeg, shapfn, shapfg, param, time, fc_q, f2e, nn, nme, nmf);
udg(:,ncu+1:end,:) = q;

Ru = cuda_Ru_elem(Ru, udg, xdg, sdg, tempen, tempeg, shapen, shapeg, param, time, fc_u, tdep, nn, nme);

Ru = cuda_Ru_face(Ru, udg, xdg, uinf, tempfn, tempfg, shapfn, shapfg, param, time, f2e, nn, nmf);
