function [udg, Ip, Im] = cuda_putfacenodes_nodal_based(udg,uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, cole2f2, ent2ind2, npf, npe, nc, e1, e2, opts)
    %TODO: STEPS
    % 1. Probably write the element-based stamping into Matlab and check...
    % 2. Jacobian stamping
    %   - We can use the "i" (DOF index or (e,n) = element+node). How to get "j"?
    %   - Probably we can start by just looping over all degrees of freedom. But then, how will we know that j is a degree of freedom that interacts with i? 
    %       - Can start with diagonal interactions where e and f are the same.
    %       - I get this intuitively but am unsure how indexing should be done. 

    % ((e,n), (e',n')) -> (f,m,l,R,L)
    % with e_R, e_L, n_R, n_L
    % n->l n'->m, R, L depend on face that i and j whether they are on the same+
    % face...sum over F is the same as the residual 
    % for residual we have (i=(e,n)) -> f, l, R, L
    
% e2 = 1082; e1 = 0;
ne = e2-e1; % number of elements
K = npf*nc; % number of degrees of freedom per face
M = npe*nc; % number of degrees of freedom per element
N = M*ne;   % total degrees of freedom in block of elements 
I = M*e1;   % starting index for residual
Ip = []; Im = [];
if opts==0
    for idx = 0:N-1
        j = rem(idx,M);               % [1, npe*nc] local dof on element
        k = (idx-j)/M+e1;               % [1, ne]     likely element
        l = rem(j,npe);               % [1, npe]    n (node on element)
        m = (j-l)/npe;                  % [1, nc]     component of current dof

        i = ent2ind1(l+npe*k+1);           %TODO: will this be off by 1? 
        if i>0
            e = i;
        else
            e = 0;
        end
%         e = (i > 0) ? i : 0;             %TODO: what is this? determine if boundary maybe? 
        n = rowe2f1((i+1)+1) - rowe2f1(e+1);   %TODO: what is this? determine # of adjacent faces? 
        for j = 0:n-1      % I think that this is looping through the faces that are adjacent to this node with positive normal vec? 
            q = cole2f1(rowe2f1(i+1)+j+1);
            p = rem(q,npf);              % [1, npf] l  (corresponding node on that face)
            s = (q-p)/npf;               % [1, nf]  f  (corresponding face on element)
            udg(I+idx+1) = udg(I+idx+1) - uh(p+npf*m+K*s+1); 
%             disp("what?")
            Im = [Im; I+idx+1, p+npf*m+K*s+1];
        end
        i = ent2ind2(l+npe*k+1);           %TODO: will this be off by 1? 
        if i>0
            e = i;
        else
            e = 0;
        end
%         e = (i > 0) ? i : 0;            %TODO: what should this be? Certainly not zero...
        n = rowe2f1((i+1)+1) - rowe2f1(e+1); 
        for j = 0:n-1  % and this will be looping through faces with negative normal perhaps? 
            q = cole2f2(rowe2f2(i+1)+j+1);
            p = rem(q,npf);              % [1, npf]
            s = (q-p)/npf;                 % [1, nf]
            udg(I+idx+1) = udg(I+idx+1) + uh(p+npf*m+K*s+1); 
            Ip = [Ip; I+idx+1, p+npf*m+K*s+1];
        end                        
    end
else
    for idx = 1:N
        j = rem(idx,M)+1;     % [1, npe*nc] local dof on element
        k = (idx-j)/M+e1;     % [1, ne]     likely element
        l = rem(j,npe)+1;     % [1, npe]    n (node on element)
        m = (j-l)/npe;        % [1, nc]     component of current dof

        i = ent2ind1(l+npe*k)+1;           %TODO: will this be off by 1? 
        if i>1
            e = i+1;
        else
            e = 1;
        end            %TODO: what is this? determine if boundary maybe? 
        n = rowe2f1(i+1) - rowe2f1(e);   %TODO: what is this? determine # of adjacent faces? 
        for j = 1:n           % I THINK that this is looping through the faces that are adjacent to this node with positive normal vec? 
            q = cole2f1(rowe2f1(i)+j);
            p = rem(q,npf)+1;              % [1, npf] l  (corresponding node on that face)
            s = (q-p)/npf;               % [1, nf]  f  (corresponding face on element)
            udg(I+idx) = udg(I+idx) - uh(p+npf*m+K*s); 
        end
    end    
end

%     a = 1;
%     nf = f2-f1+1;
%     ndf = npf*nf;

    %NOTE: This is based on the face-based assembly...we should move to the element-based form used by Exasim

    
    %% ELEMENT-BASED PUTFACENODES
% template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
%         int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts)
% {    
%     int ne = e2-e1;
%     int K = npf*nc;
%     int M = npe*nc;
%     int N = M*ne;        
%     int I = M*e1; // starting point of indexing; remember we're going through blocks of elements for this version, not faces
%     if (opts==0) {
%         for (int idx = 0; idx<N; idx++) { // go through all degrees of freedom
%             int j = idx%M;              //[1, npe*nc] 
%             int k = (idx-j)/M+e1;       //[1, ne]     // e (element) 
%             int l = j%npe;              //[1, npe]    // n (node on element)
%             int m = (j-l)/npe;          //[1, nc]     // ncu
%             int q, p, s;
            
%             int i = ent2ind1[l+npe*k];
%             int e = (i > 0) ? i : 0;
%             int n = rowe2f1[i+1] - rowe2f1[e];
%             for (j=0; j<n; j++) {  // // I THINK that this is looping through the faces that are adjacent to this node with positive normal vec? 
%                 q = cole2f1[rowe2f1[i]+j];
%                 p = q%npf;              // [1, npf] l  (corresponding node on that face)
%                 s = (q-p)/npf;          // [1, nf]  f  (corresponding face on element)
%                 udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
%             }            
            
%             i = ent2ind2[l+npe*k];
%             e = (i > 0) ? i : 0;
%             n = rowe2f2[i+1] - rowe2f2[e];
%             for (j=0; j<n; j++) {  // and this will be looping through
%                 q = cole2f2[rowe2f2[i]+j];
%                 p = q%npf;              // [1, npf]
%                 s = (q-p)/npf;          // [1, nf]      
%                 udg[I+idx] = udg[I+idx] + uh[p+npf*m+K*s]; 
%             }                        
%         }        
%     }
%     else {
%         for (int idx = 0; idx<N; idx++) {
%             int j = idx%M;              //[1, npe*nc]
%             int k = (idx-j)/M+e1;       //[1, ne]      
%             int l = j%npe;              //[1, npe]
%             int m = (j-l)/npe;          //[1, nc] 
            
%             int i = ent2ind1[l+npe*k];
%             int e = (i > 0) ? i : 0;
%             int n = rowe2f1[i+1] - rowe2f1[e];
%             for (j=0; j<n; j++) {
%                 int q = cole2f1[rowe2f1[i]+j];
%                 int p = q%npf;          // [1, npf]
%                 int s = (q-p)/npf;          // [1, nf]    
%                 udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
%             }            
%         }
%     }
% }

    