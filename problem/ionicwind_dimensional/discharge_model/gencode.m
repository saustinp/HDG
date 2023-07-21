function gencode(p, f, udg, param, time, fname_prefix, ndim, nc_fhat)
syms zero one

filename1 = [fname_prefix num2str(ndim) 'd.m'];

%%% compute Jacobian
jac_f  = jacobian(f,udg);

%%% And patch with vector zero or one to have the right sizes
for ii = 1:size(f,1)
    for jj = 1:size(f,2)
        temp = ismember(symvar(f(ii,jj)),udg);
        if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
        elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
        end
    end
end

for ii = 1:size(jac_f,1)
    for jj = 1:size(jac_f,2)
        temp = ismember(symvar(jac_f(ii,jj)),udg);
        if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
        elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
        end
    end
end    

% generate a temporary matlab file      % NOTE: "p" was changed from the template
matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{p,udg,param,time,[zero one]},'outputs', {'out','jac_out'});

% open the file and modify it
fid = fopen('tmp.m','r');
gid = fopen(filename1,'wt');

tline = fgetl(fid); 
i=1;       
end_count = 0;
while ischar(tline)        
    str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
    str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
    str = strrep(str, 'in1', 'p');          % NOTE: This was changed from the template
    str = strrep(str, 'IN1', 'P');        
    str = strrep(str, 'in2', 'udg');
    str = strrep(str, 'IN2', 'UDG');        
    str = strrep(str, 'in3', 'param');
    str = strrep(str, 'IN3', 'PARAM');                
    str = strrep(str, ',in5)', ')');                
    str = strrep(str, ',IN5)', ')');    
    str = strrep(str, 'param(:,1)', 'param{1}'); 
    str = strrep(str, 'param(:,2)', 'param{2}'); 
    str = strrep(str, 'param(:,3)', 'param{3}'); 
    str = strrep(str, 'param(:,4)', 'param{4}'); 
    str = strrep(str, 'param(:,5)', 'param{5}'); 
    str = strrep(str, 'param(:,6)', 'param{6}'); 
    str = strrep(str, 'param(:,7)', 'param{7}'); 
    str = strrep(str, 'param(:,8)', 'param{8}'); 
    str = strrep(str, 'param(:,9)', 'param{9}');     
    str = strrep(str, 'param(:,10)', 'param{10}');      % NOTE: These lines were added
    str = strrep(str, 'param(:,11)', 'param{11}');     
    str = strrep(str, 'param(:,12)', 'param{12}');     
    str = strrep(str, 'param(:,13)', 'param{13}');     
    str = strrep(str, 'param(:,14)', 'param{14}');     
    str = strrep(str, 'param(:,15)', 'param{15}');     
    str = strrep(str, 'param(:,16)', 'param{16}');     
    str = strrep(str, 'param(:,17)', 'param{17}');     
    str = strrep(str, 'param(:,18)', 'param{18}');     
    str = strrep(str, 'param(:,19)', 'param{19}');     
    str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
    str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
    if i==7
        str = '[ng,nc] = size(udg);';
        fprintf(gid, '%s\n', str);                  
        str = ['nch = ' num2str(nc_fhat) ';'];
        fprintf(gid, '%s\n', str);                  
        str = ['nd = ' num2str(ndim) ';'];
    end

    % Strips last 'end' form function and moves it to the very end
    if contains(str, 'end')
        end_count = end_count + 1;
    end
    if end_count == 2
        str = strrep(str, 'end', '');
    end
    fprintf(gid, '%s\n', str);                  
    tline = fgetl(fid);        
    i=i+1;
end

% NOTE: This will produce a function that has errors. You need to move the following "reshape" lies above the second "end" statement in order to be enclosed in the the function
if fname_prefix == "flux"
    fprintf(gid, '%s\n', 'out = reshape(out,ng,nch,nd);');
    fprintf(gid, '%s\n', 'jac_out = reshape(jac_out,ng,nch,nd,nc);'); 
elseif fname_prefix == "source"
    fprintf(gid, '%s\n', 'out = reshape(out,ng,nch);');
    fprintf(gid, '%s\n', 'jac_out = reshape(jac_out,ng,nch,nc);');      
end            

fprintf(gid, '%s\n', 'end');                  

fclose(fid);
fclose(gid);
delete('tmp.m');