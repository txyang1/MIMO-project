function [Grad] = MSEDualMACgrad(H,Qblk)
% function [Grad] = MSEDualMACgrad(H,Qblk)
%
% MSEDualMACgrad computes the gradient evaluated at Qblk. The result is
% returned as a block diagonal matrix with the gradients in each block
% Input:
% H: K x 1 cell array with channels H_k
% Qblk: KM x KM block diagonal matrix with transmit covariance matrices
% Output:
% Grad: KM x KM block diagonal matrix with gradient in each block

K = length(H);
M = size(H{1},1);
N = size(H{1},2);

I  = eye(N);

Q = mat2cell(Qblk,[M M],[M M]);

X = zeros(N,N);
for i = 1:K
    X = X + H{i}'*Q{i,i}*H{i};
end
invX = inv(I+X);

tempgrad = cell(K,1);
for l = 1:K
    tempgrad{l}  = -invX*H{l}'*H{l}*invX;
end
Grad = blkdiag(tempgrad{:});

end

