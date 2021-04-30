function [L,S,obj,err,iter] = trpca_tnn(X,lambda,transform,opts)

% Tensor robust principal component analysis based on tensor nuclear norm 
% under linear transform
%
% min_{L,S} ||L||_*+lambda*||S||_1, s.t. X=L+S
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       T       -    d3*d3 matrix, T'*T=T*T'=l*I, l>0
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%   transform   - a structure which defines the linear transform
%       transform.L: the linear transform of two types:
%                  - type I: function handle, i.e., @fft, @dct
%                  - type II: invertible matrix of size n3*n3
%
%       transform.inverseL: the inverse linear transform of transform.L
%                         - type I: function handle, i.e., @ifft, @idct
%                         - type II: inverse matrix of transform.L
%
%       transform.l: a constant which indicates whether the following property holds for the linear transform or not:
%                    L'*L=L*L'=l*I, for some l>0.                           
%                  - transform.l > 0: indicates that the above property holds. Then we set transform.l = l.
%                  - transform.l < 0: indicates that the above property does not hold. Then we can set transform.l = c, for any constant c < 0.
%
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations
%
% version 1.0 - 04/03/2019
% version 1.1 - 29/04/2021
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
% References:
% Canyi Lu, Pan Zhou, Exact recovery of tensor robust principal component 
% analysis under linear transforms. arXiv preprint arXiv:1907.08288, 2019
%

tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if nargin < 3
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

dim = size(X);
L = zeros(dim);
S = L;
Y = L;

for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % update L
    [L,tnnL] = prox_tnn(-S+X-Y/mu,1/mu,transform);
    % L = real(L); % uncomment this line if transform dftmtx(n3) is used to avoid very small complex element

    % update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);

    dY = L+S-X;
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnL+lambda*norm(S(:),1);
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err) ]);
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);
end
obj = tnnL+lambda*norm(S(:),1);
err = norm(dY(:));
