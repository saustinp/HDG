classdef mArray  
    methods
%         function array = mArray(sizes,indices)
%             %array.theSizes = sizes;
%             %array.theIndices = indices;
%         end

%         function array = mArray()
%             array.theArray = [];
%         end
               
        function array = mArray(A)
            if nargin > 0
                if isa(A,'mArray')
                    array.theArray = A.theArray;
                else
                    array.theArray = A;
                end
            else
                array.theArray = [];
            end
        end
        
        function indices = getIndices(A)
            indices = A.theIndices(A.theIndices ~= 1);
        end
        
        function remainingIndices = substractIndices(A,indices)
            allIndices = getIndices(A);
            remainingIndices = allIndices(allIndices ~= indices);
        end
        
        function S2 = getOrderedSubs(A,S)
            nIndices = numel(S.subs);
            subs = S.subs;
            
            nAllIndices = ndims(A.theArray);
            
            S2.type = '()';
            S2.subs = {};
            
            for k = 1:nAllIndices
                S2.subs{k} = ':';
            end
            
            for k = 1:nIndices / 2
                pos = subs{1,2*k -1};
                range = subs{1,2*k};
                S2.subs{pos} = range;
            end
            
        end
        
        function B = subsref(A,S)
            switch S.type
                case '.'
                    A = builtin('subsref',A,S);
                case '()'
%                     A = mArray(A);
                    S2 = getOrderedSubs(A,S);
                    B = mArray();
                    B.theArray = subsref(A.theArray,S2);
            end
        end
        
        function A = subsasgn(A,S,B)
            switch S.type
                case '.'
                    A = builtin('subsasgn',A,S,B);
                case '()'
                    S2 = getOrderedSubs(A,S);
%                     A = mArray(A);
                    B = mArray(B);
                    A.theArray = subsasgn(A.theArray,S2,B.theArray); % It is required
            end
        end
        
        function [] = disp(A)
            disp(A.theArray);
        end
        
        function value = size(A,varargin) % TODO put the additional parameters
            if(nargin > 1)
                value = size(A.theArray,varargin{:});
            else
                value = size(A.theArray);
            end
        end
        
        function value = numel(A)
            value = numel(A.theArray);
        end
        
        function text = indices(A,labels)
            sizeA = size(A);
            sizes = sizeA(size(A) ~= 1);
            value = labels(size(A) ~= 1);
            text = {};            
            text(1,:) = value;
            text(2,:) = num2cell(sizes); 
        end
        
        function [value] = sizes(A)
            sA = size(A);
            value = sA(sA > 1);
        end
        
        function [nDims] = ndims(A)
            nDims = ndims(A.theArray);
        end
        
        function [data] = getData(A)
            data = A.theArray;
        end
        
        function [C] = cat(dim,varargin)            
            nArrays = nargin - 1;
            arrays = cell(1,nArrays);            
            for iArray = 1:nArrays
                array = varargin{iArray};
                arrays{iArray} = array.theArray;
            end
            data = cat(dim,arrays{:});
            C = mArray(data);
        end
        
        function B = map(A,ii,jj)
%             ii = sort(ii);
%             jj = sort(jj);
            maxIndex = max([max(ii) max(jj) max(ndims(A.theArray))]);
            permutation = 1:maxIndex;
            
            permutation(jj) = ii;
            permutation(ii) = jj;
            
            B = mArray();
            
            B.theArray = permute(A.theArray,permutation);
        end
        
        function B = join(A,indices,index)
%             indices = sort(indices);
            sizesA = size(A);
            maxIndex = max(max([indices index]),numel(sizesA));
            nIndices = numel(indices);
            
            allIndices = 1:maxIndex;
            
            indicesToOffset = allIndices(allIndices >= index);
                       
            offset = zeros(1,maxIndex);                       
            offset(indicesToOffset) = nIndices;
                                               
            ii = offset(indices) + indices;
            jj = index:index+(nIndices-1);
                       
            sizesD = ones(1,maxIndex + nIndices);
            sizesD(allIndices + offset) = sizesA(allIndices);
            sizesD(jj) = 1;
            
            Ddata = reshape(A.theArray,sizesD);
            
            D = mArray(Ddata);            
            
            C = map(D,ii,jj);
                        
            sizesC = sizesA;
            
            sizesC(indices) = 1;
            sizesC(index) = prod(sizesA(indices));
            
            data = reshape(C.theArray,sizesC);            
            B = mArray(data);              
        end
        
        function B = split(A,index,indices,sizes)
            %indices = sort(indices);            
            maxIndex = max(max([index indices]),ndims(A.theArray));
            nIndices = numel(indices);
            
            allIndices = 1:maxIndex;
            
            sizesA = ones(1,maxIndex+nIndices);
            sizesA(1:ndims(A.theArray)) = size(A);
            
            indicesToOffset = allIndices(allIndices >= index);
                       
            offset = zeros(1,maxIndex);                       
            offset(indicesToOffset) = nIndices;
               
            ii = index:index+(nIndices-1);
            jj = offset(indices) + indices;
            
            sizesD = ones(1,maxIndex + nIndices);
            sizesD(allIndices + offset) = sizesA(allIndices);
            sizesD(ii) = sizes;
            sizesD(index+nIndices) = 1;
                        
            Ddata = reshape(A.theArray,sizesD);
            
            D = mArray(Ddata);             
            
            C = map(D,ii,jj); 
            
            sizesC = ones(1,maxIndex+nIndices);
            sizesC(1:numel(size(C))) = size(C);

            sizesB = sizesC([1:index-1 (index+nIndices:maxIndex+nIndices)]);
                        
            data = reshape(C.theArray,sizesB);            
            B = mArray(data);
            
        end
        
        function [C] = sum(A,k)
            
            if(isempty(k))
                C = A;
                return
            end
            
            if ~isa(A,'mArray')
                A = mArray(A);
            end
            
            k = sort(k);
            
            sA = size(A);
            maxIndex = max([numel(sA),max(k)]);
            ind = 1:maxIndex; 
            
            sA = ones(1,maxIndex);
            
            sizeA = size(A);
            
            sA(1:numel(sizeA)) = size(A);
            
            rA = sA;
            
            rA(k) = 0;
            
            l = ind( rA >  1);
            z = ind( rA == 1);
                                    
            q = sA(k);
            r = sA(l);
            s = sA(z);
            
            qq = prod(q);
            rr = prod(r);
            ss = prod(s);
            
            Ad = A.theArray;
            Ap = permute(Ad,[k l z]);
            Ar = reshape(Ap,[qq rr ss 1 1]);
            Cc = sum(Ar,1);
            Cr = reshape(Cc,[ones(1,numel(k)) r s 1 1]);
            Cp = ipermute(Cr,[k l z]);
            C = mArray(Cp);
        end
                
        function [C,ops,mops] = contract(A,B,k)
            if ~isa(A,'mArray')
                A = mArray(A);
            end
            
            if ~isa(B,'mArray')
                B = mArray(B);
            end
            
            k = sort(k); % Sort k allows to reduce unrequired permutations
            
            sA = size(A);
            sB = size(B);
            maxIndex = max([numel(sA),numel(sB),max(k)]);
            ind = 1:maxIndex; 
            
            sA = ones(1,maxIndex);
            sB = ones(1,maxIndex);
            
            sizeA = size(A);
            sizeB = size(B);
            
            sA(1:numel(sizeA)) = size(A);
            sB(1:numel(sizeB)) = size(B);
            
            rA = sA;
            rB = sB;
            
            rA(k) = 0;
            rB(k) = 0;
            
            ki = ind(sA > 1 & sA ~= sB & rA == 0); % Detect non-singleton indices to contract in i that are not in j
            kj = ind(sB > 1 & sA ~= sB & rB == 0); % Detect non-singleton indices to contract in j that are not in i
            
            rA(ki) = 1;
            rB(kj) = 1;
                      
            k = ind(rA == 0 & rB == 0); % Put in k just those indices that have to be contracted both in A and B
            
            rA(kj) = 1;
            rB(ki) = 1;
            
            l = ind(rA == rB & rA >  1 & rB >  1);
            z = ind(rA == rB & rA == 1 & rB == 1); % Detect singletons and ki and kj (because the are going to be singletons before contractBase, allows to remove singletos avoiding unrequired permutations
            
            rA(l) = 0;
            rB(l) = 0;
            
            i = ind(rA > 1);
            j = ind(rB > 1);
            
            m = sA(i);
            n = sB(j);
            q = sA(k);
            r = sA(l);
            s = sA(z);
            
            mm = prod(m);
            nn = prod(n);
            qq = prod(q);
            rr = prod(r);
            ss = prod(s);
            
            As = sum(A,ki); % Contract ki, allows to preserve contraction code when A does not depend on kj
            Bs = sum(B,kj); % Contract kj, allows to preserve contraction code when B does not depend on ki
                                                                                  
            Cr = contractBase(As,Bs,i,j,l,k,z,m,n,r,s,mm,nn,qq,rr,ss);
            
            C = mArray(Cr);
                        
            ops = 2*mm*nn*qq*rr;
            mops = numel(C) + numel(A) + numel(B);
        end
        
        function Cs = contractBase(A,B,i,j,l,k,z,m,n,r,s,mm,nn,qq,rr,ss)
            % \warning May be I sould sort i j l to improve performance
            maxIndex = max([i k j l z]);
            isz = zeros(1,maxIndex);
            isk = zeros(1,maxIndex);
            isz(z) = 1;
            isk(k) = 1;
            
            squeezeOffset = -cumsum(isz);  % Shift offset indices due to singletons
            squeezeOffsetK = -cumsum(isk); % Shift offset indices due to contraction of k
                        
            % Shifting indices due to singletons
            is = i + squeezeOffset(i);
            js = j + squeezeOffset(j);
            ks = k + squeezeOffset(k);
            ls = l + squeezeOffset(l);
            
            % Shifting indices due to contractio of k, for final result C
            it = is + squeezeOffsetK(i);
            jt = js + squeezeOffsetK(j);
            lt = ls + squeezeOffsetK(l);
            kt = ks + squeezeOffsetK(k);
            
            sA = ones(1,maxIndex);
            sB = ones(1,maxIndex);
            
            sA(1:numel(size(A))) = size(A);
            sB(1:numel(size(B))) = size(B);
            
            sC = max(sA,sB);
            sC(k) = 1; % Preparing sizes for final result. It does not depend on contain k
            
            % For removing singletons in A and B
            rA = sA(isz == 0); % For removing singletons in A
            rB = sB(isz == 0); % For removing singletons in B
            
            % Removing singletons in A and B by means of rA and rB
            Ad = reshape(A.theArray,[rA 1 1]); 
            Bd = reshape(B.theArray,[rB 1 1]); 
            
                Ap = mPermute(Ad,[is ks js ls]);
                Bp = mPermute(Bd,[is ks js ls]);             

                As = reshape(Ap,[mm qq rr 1 1]);
                Bs = reshape(Bp,[qq nn rr 1 1]);
                
                Cs = zeros([mm nn rr]);
                
                for ll = 1:rr
                    Cs(:,:,ll) = As(:,:,ll) * Bs(:,:,ll);
                end
                               
                Cp = reshape(Cs, [m  n  r 1 1]);
                order = [it jt lt];
                if(numel(order) > 1)
                    Ct = ipermute(Cp,[it jt lt   ]);
                else
                    Ct = Cp;
                end
                Cr = reshape(Ct,sC);
                
                Cs = reshape(Cr,sC); % Reshape to proper size     
                
                return
            
            if(mm == 2 && qq == 2 && nn == 2) % These specializations could be moved               
                Ap = permute(Ad,[ls is ks js]);
                Bp = permute(Bd,[ls is ks js]);

                As = reshape(Ap,[rr mm qq 1 1]);
                Bs = reshape(Bp,[rr qq nn 1 1]); 
                
                Cs = zeros([rr mm nn]);               % to functions and decide on contract which one you call
                Cs(:,1,1) = As(:,1,1) .* Bs(:,1,1) + As(:,1,2) .* Bs(:,2,1);
                Cs(:,1,2) = As(:,1,1) .* Bs(:,1,2) + As(:,1,2) .* Bs(:,2,2);
                Cs(:,2,1) = As(:,2,1) .* Bs(:,1,1) + As(:,2,2) .* Bs(:,2,1);
                Cs(:,2,2) = As(:,2,1) .* Bs(:,1,2) + As(:,2,2) .* Bs(:,2,2);
                
                Cp = reshape(Cs, [r  m  n 1 1]);
                Ct = ipermute(Cp,[lt it jt   ]);
                Cr = reshape(Ct,sC);
            elseif (mm < 32 && qq < 32 && nn < 32 && rr > 1)
                Ap = permute(Ad,[ls is ks js]);
                Bp = permute(Bd,[ls is ks js]);             

                As = reshape(Ap,[rr mm qq 1 1]);
                Bs = reshape(Bp,[rr qq nn 1 1]);
                
                Cs = zeros([rr mm nn]);
                        
                for jj = 1:nn
                    for ii = 1:mm
                        c = zeros(rr,1);
                        for kk = 1:qq
                            c = c + As(:,ii,kk) .* Bs(:,kk,jj);
                        end
                        Cs(:,ii,jj) = c;
                    end
                end
                               
                Cp = reshape(Cs, [r m  n  1 1]);
                order = [lt it jt];
                if(numel(order) > 1) % Check for avoiding permutation of only 1 index
                    Ct = ipermute(Cp,[lt it jt]);
                else
                    Ct = Cp;
                end
                Cr = reshape(Ct,sC);
            elseif (rr==1)
                if(numel([is ks js ls]) > 1) % Check for avoiding permutation of only 1 index
                    Ap = permute(Ad,[is ks js ls]);
                    Bp = permute(Bd,[is ks js ls]);
                else
                    Ap = Ad;
                    Bp = Bd;
                end
                                
                As = reshape(Ap,[mm qq rr 1 1]);
                Bs = reshape(Bp,[qq nn rr 1 1]);
                
                Cs = As * Bs;
      
                Cp = reshape(Cs, [m  n  r 1 1]);
                %Ct = ipermute(Cp,[it jt lt   ]);
                order = [it jt lt];
                if(numel(order) > 1)
                    Ct = ipermute(Cp,order);
                else
                    Ct = Cp;
                end
                Cr = reshape(Ct,sC);           
            else
                Ap = permute(Ad,[is ks js ls]);
                Bp = permute(Bd,[is ks js ls]);             

                As = reshape(Ap,[mm qq rr 1 1]);
                Bs = reshape(Bp,[qq nn rr 1 1]);
                
                Cs = zeros([mm nn rr]);
                for ll = 1:rr
                    Cs(:,:,ll) = As(:,:,ll) * Bs(:,:,ll);
                end
                               
                Cp = reshape(Cs, [m  n  r 1 1]);
                order = [it jt lt];
                if(numel(order) > 1)
                    Ct = ipermute(Cp,[it jt lt   ]);
                else
                    Ct = Cp;
                end
                Cr = reshape(Ct,sC);
            end

                   
            Cs = reshape(Cr,sC); % Reshape to proper size          
        end
        
        function Cr = contractBase3(A,B,i,j,l,k,m,n,r,mm,nn,qq,rr)
            Ap = permute(A.theArray,[i k l j]);
            Bp = permute(B.theArray,[k j l i]);
            
            As = reshape(Ap,[mm qq rr 1]);
            Bs = reshape(Bp,[qq nn rr 1]);            
            
%             Cs = zeros([mm nn rr]);
            
            Cs = mContract(As,Bs);
            
            Cp = reshape(Cs,[m n r 1 1]);
            
            Cr = ipermute(Cp,[i j l k]);
        end
        
        function Cr = contractBase2(A,B,i,j,l,k,m,n,r,mm,nn,qq,rr)
            Ap = permute(A.theArray,[l i k j]);
            Bp = permute(B.theArray,[l j k i]);
            
            As = reshape(Ap,[rr mm 1  qq 1 1]);
            Bs = reshape(Bp,[rr 1  nn qq 1 1]);            
            
            Cs = zeros([rr mm nn]);
            
%             for ll = 1:rr
%                 Cs(:,:,ll) = As(:,:,ll) * Bs(:,:,ll);
%             end

%             for jj = 1:nn
%                 for kk = 1:qq
%                     Bsjjkk = Bs(:,1,jj,kk);
%                     for ii = 1:mm                        
%                         Cs(:,ii,jj) = Cs(:,ii,jj) + As(:,ii,1,kk) .* Bsjjkk;                
%                     end
%                 end
%             end

            for ii = 1:mm
                for jj = 1:nn
                    Csiijj = zeros(rr,1);
                    for kk = 1:qq
                        Csiijj = Csiijj + As(:,ii,1,kk) .* Bs(:,1,jj,kk);
                    end
                    
                    Cs(:,ii,jj) = Csiijj;
                end
            end
            
            Cp = reshape(Cs,[m n r 1 1]);
            
            Cr = ipermute(Cp,[i j l k]);            
        end
        
        function C = mFun(func,A,B)
            A = mArray(A);
            B = mArray(B);
            C = mArray();
            C.theArray = bsxfun(func,A.theArray,B.theArray);
        end
        
        function C = plus(A,B)
            C = mFun(@plus,A,B);
        end
        
        function C = minus(A,B)
            C = mFun(@minus,A,B);

        end
        
        function C = uminus(A)
            A = mArray(A);
            C = mArray();
            C.theArray = uminus(A.theArray);
        end
        
        function C = uplus(A)
            A = mArray(A);
            C = mArray();
            C.theArray = uplus(A.theArray);
        end
        
        function C = times(A,B)
            C = mFun(@times,A,B);           
        end
        
        function C = mtimes(A,B)
            C = mFun(@times,A,B);  
        end
        
        % To do mtimes (maybe to contract repeated indices)
        
        function C = rdivide(A,B)
            C = mFun(@rdivide,A,B);  
        end
        
        function C = ldivide(A,B)
            C = mFun(@ldivide,A,B);
        end
        
        function C = power(A,B)
            C = mFun(@power,A,B);            
        end
        
        function C = lt(A,B)
            C = mFun(@lt,A,B); 
        end
        
        function C = gt(A,B)
            C = mFun(@gt,A,B); 
        end
        
        function C = le(A,B)
            C = mFun(@le,A,B);
        end        
                
        function C = ge(A,B)
            C = mFun(@ge,A,B);
        end  
                
        function C = ne(A,B)
            C = mFun(@ne,A,B);
        end  
                
        function C = eq(A,B)
            C = mFun(@eq,A,B);
        end          
                
        function C = and(A,B)
            C = mFun(@and,A,B);
        end          
                
        function C = or(A,B)
            C = mFun(@or,A,B);
        end         
                
        function C = not(A,B)
            C = mFun(@not,A,B);
        end
        
        function B = double(A)
            B = mArray(A);
%             B.theArray = A;
        end
    end
    
    methods
        function B = sin(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = sin(A.theArray);
        end
        
        function B = cos(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = cos(A.theArray);
        end
        
        function B = tan(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = tan(A.theArray);
        end
        
        function B = exp(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = exp(A.theArray);
        end
        
        function B = sqrt(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = sqrt(A.theArray);
        end
        
        function B = abs(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = abs(A.theArray);
        end
        
        function B = sign(A)
            % A = mArray(A);
            B = mArray();            
            B.theArray = sign(A.theArray);
        end
%        function B = sum(A,indices)
%            % A = mArray(A);
%            B = mArray(A);
%            nIndices = numel(indices);
%            for iIndex = 1:nIndices
%                B.theArray = sum(B.theArray,indices(iIndex));
%            end
%        end
    end
    
    methods(Static)
        
        function [tokens] = tokenizeLabels(labels)
            tokens = {};
            remain = labels;
            while (~isempty(remain))
                [token,remain] = strtok(remain);
                tokens = [tokens token];
            end            
        end

        function [p] = initIndices(positions,sizes,labels,captions)
            nIndices = numel(positions);

            tokens = mArray.tokenizeLabels(labels);

            for iIndex = 1:nIndices
                pString = ['p.'  tokens{iIndex} '=' int2str(positions(iIndex)) ';'];
                eval(pString);
            end
            
            for iIndex = 1:nIndices
                nString = ['p.n' tokens{iIndex} '=' int2str(sizes(iIndex)) ';'];
                eval(nString);
            end

            p.n = sizes;
            p.labels = tokens;
            p.positions = positions;
            if (exist('captions','var'))
                p.captions = captions;
            end

            % p
            % p.n
            % p.labels{:}
            % p.captions{:}
        end
        
        function array = zeros(sizes,indices)
            array = mArray(); % mArray(sizes,indices);
            array.theArray = zeros(mArray.getSizes(sizes,indices));
        end
        
        function array = ones(sizes,indices)
            array = mArray(); % mArray(sizes,indices);
            array.theArray = ones(mArray.getSizes(sizes,indices));
        end
        
        function array = rand(sizes,indices)
            array = mArray(); % mArray(sizes,indices);
            array.theArray = rand(mArray.getSizes(sizes,indices));
        end
        
        function s = getSizes(sizes,indices)
            maxIndex = max(indices)+1;
            s = ones(1,maxIndex);
            s(indices) = sizes;
        end
        
        function array = fromArray(A,indices) % Indices come in order of A
            Ap = inversePermute(A,indices,max(indices));
            
            array = mArray(Ap);
        end
        
        function B = toArray(A,indices)            
            B = mPermute(A.theArray,indices);            
        end
        
        function B = zerosTransform(sizes,reducedIndices,expandedIndices)
            B = mArray();
            n = sizes;
            n(reducedIndices) = 1;
            n(expandedIndices) = sizes(expandedIndices);
            
            B.theArray = zeros(n);
        end
        
    end
    
    properties
        %theSizes;
        %theIndices;
        theArray;
    end
    
end

