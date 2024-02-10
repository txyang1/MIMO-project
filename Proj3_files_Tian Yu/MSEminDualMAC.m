function [W,G] = MSEminDualMAC(H, Ptx)
% function [W,G] = MSEminDualMAC(H, Ptx)
%
% The function computes the MSE optimal transmit and receive filters the
% dual multiple access channel H subject to a sum power constraint Ptx. A
% gradient projection algorithm is used to compute the solution
% Input:
% H: K x 1 cell array with user channels
% Ptx: Available sum transmit power
% Output:
% W: K x 1 cell array with transmit filters
% G: K x 1 cell array with receive filters


K = length(H);
[M,N] = size(H{1});

% Initialize block diagonal transmit covariance
QblkInit = eye(M*K)*Ptx/(K*M);

% Define function handels for the gradient projection algorithm
fun = @(X) MSEDualMAC(H,X);
grad = @(X) MSEDualMACgrad(H,X);
proj = @(X) MSEDualMACproj(X,Ptx);

% Run gradient projection algorithm
[Qopt,~,~] = projGrad(fun, grad, proj, QblkInit, 100);

% Compute optimal filters from block diagonal transmit covariance Qopt
Hco = vertcat(H{:})';
X = eye(N) + Hco*Qopt*Hco';

Q = cell(K,1);
W = cell(K,1);
G = cell(K,1);
for k=1:K
    Q{k} = Qopt(M*(k-1)+1:M*k,M*(k-1)+1:M*k);
    [U,S,~] = svd(Q{k});
    W{k} = U*sqrt(S);
    G{k} = W{k}'*H{k}*inv(X);
end

end

