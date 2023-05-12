function udg = cuda_putfacenodes_Jac(udg,up,um,f2e,npf,ncu,npe,f1,f2,opts)

    a = 1;
    nf = f2-f1+1;
    ndf = npf*nf;
    %TODO: STEPS
    % 1. Probably write the element-based stamping into Matlab and check...
    % 2. Jacobian stamping
    %   - We can use the "i" (DOF index or (e,n) = element+node). How to get "j"?
    %   - Probably we can start by just looping over all degrees of freedom. But then, how will we know that j is a degree of freedom that interacts with i? 
    %       - Can start with diagonal interactions where e and f are the same.
    %       - I get this intuitively but am unsure how indexing should be done. 

    % ((e,n), (e',n')) -> (f,m,l,R,L)
    % with e_R, e_L, n_R, n_L
    % n->l n'->m, R, L depend on face that i and j whether they are on the same
    % face...sum over F is the same as the residual 
    % for residual we have (i=(e,n)) -> f, l, R, L
    
    %NOTE: This is based on the face-based assembly...we should move to the element-based form used by Exasim
    udg = reshape(udg, [ncu*npe*1082, ncu*npe*1082]);
    if opts==0
        for i=1:ndf
            m = npf*(f1-1)+i;
            k1 = f2e(1,m); m1 = rem(k1-1,npe)+1; n1 = (k1-m1)/npe+1;
            k2 = f2e(2,m); m2 = rem(k2-1,npe)+1; n2 = (k2-m2)/npe+1;      
            for j=1:ncu
                Ik1 = m1+(j-1)*npe+(n1-1)*npe*ncu;
                Ik2 = m2+(j-1)*npe+(n2-1)*npe*ncu;
                %TODO: this must be wrong...just curious about indexing....
                %TODO: how should up and um be arranged and indexded into?
                udg(Ik1, Ik1) = udg(Ik1, Ik1) - a*up(i+(j-1)*ndf); 
                udg(Ik1, Ik2) = udg(Ik1, Ik2) - a*up(i+(j-1)*ndf); 
    
                udg(Ik2, Ik1) = udg(Ik2, Ik1) + a*um(i+(j-1)*ndf);  
                udg(Ik2, Ik2) = udg(Ik2, Ik2) + a*um(i+(j-1)*ndf);
            end
        end
    else
        for i=1:ndf
            m = npf*(f1-1)+i;
            k1 = f2e(1,m); m1 = rem(k1-1,npe)+1; n1 = (k1-m1)/npe+1;       
            for j=1:ncu
                udg(m1+(j-1)*npe+(n1-1)*npe*ncu) = udg(m1+(j-1)*npe+(n1-1)*npe*ncu) - a*uh(i+(j-1)*ndf); 
            end
        end    
    end
    
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
    
    % function udg = cuda_putfacenodes(udg,uh,f2e,npf,ncu,nde,f1,f2,opts)
    % 
    % a = 1;
    % nf = f2-f1+1;
    % ndf = npf*nf;
    % if opts==0
    %     for i=1:ndf
    %         m = npf*(f1-1)+i;
    %         k1 = f2e(1,m);
    %         k2 = f2e(2,m);    
    %         for j=1:ncu
    %             udg(k1+(j-1)*nde) = udg(k1+(j-1)*nde) - a*uh(i+(j-1)*ndf); 
    %             udg(k2+(j-1)*nde) = udg(k2+(j-1)*nde) + a*uh(i+(j-1)*ndf);  
    %         end
    %     end
    % else
    %     for i=1:ndf
    %         m = npf*(f1-1)+i;
    %         k1 = f2e(1,m);        
    %         for j=1:ncu
    %             udg(k1+(j-1)*nde) = udg(k1+(j-1)*nde) - a*uh(i+(j-1)*ndf); 
    %         end
    %     end    
    % end
    % 
    