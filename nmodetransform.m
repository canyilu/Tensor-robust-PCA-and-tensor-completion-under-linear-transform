function B = nmodetransform(A,transform,n)
% Calculates the transform on the n-Mode of a Tensor A, the size should keep
% the same
% B = nmodetransform(A,fun,n)
% 
% B = A (x)_n M .. According to the Definition in De Lathauwer (2000)
%
% with:
% A:    (I_1 x I_2 x .. I_n x .. I_N) .. ->  n is in [1..N]
% M:    (J   x I_n)
% B:    (I_1 x I_2 x .. J x   .. I_N)
%
% note: "(x)_n" is the operator between the tensor and the matrix
% 
% 

% check inputs:
dimvec = size(A);
n = fix(n);

if (length(dimvec)<n || n<1)
    error('nmodeproduct: n is not within the order range of tensor A ');
end

% shift A to prepare flattening: (i.e. make dimension 1 (columns) to 'n', the one we would like to replace)
Ash = shiftdim(A,n-1);
dimvecB = size(Ash);

if strcmp(transform,'fwt')
    J = log2(dimvecB(1));
    B = fwt(Ash(:,:),'db8',J);
elseif strcmp(transform,'ifwt') == 1
    J = log2(dimvecB(1));
    B = ifwt(Ash(:,:),'db8',J,dimvecB(1));
elseif strcmp(transform,'ufwt')
    J = log2(dimvecB(1));
    B = ufwt(Ash(:,:),'db8',J);
elseif strcmp(transform,'iufwt')
    J = log2(dimvecB(1));
    B = iufwt(Ash(:,:),'db8',J,dimvecB(1));
end

% wrap the flattened vector-array back into the previously saved tensor shape
B = reshape(B,dimvecB);

% shift the dimensions back! so that only dimension n has changed from I_n to J
B = shiftdim(B,length(dimvecB)-n+1);
