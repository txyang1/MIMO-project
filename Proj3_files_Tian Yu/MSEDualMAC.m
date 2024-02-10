function [obj] = MSEDualMAC(H, Qblk)
% function [obj] = MSEDualMAC(H, Qblk)
%
% MSEDualMAC calculates the value of the objective function in the dual
% multiple access channel.
% Input:
% H: K x 1 cell array with channels H_k
% Qblk: KM x KM block diagonal matrix with transmit covariance matrices
% 
% Output:
% obj: (real) value of the objective function

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
obj = real(trace(invX));

end

