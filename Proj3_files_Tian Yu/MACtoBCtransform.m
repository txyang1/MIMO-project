function [S] = MACtoBCtransform(Q,H,order)
% Function [S] = MACtoBCtransform(Q,H,order)
% 
% Rate based duality transformation from the
% multiple access channel to the broadcast channel.
%
% Inputs
% Q: K x 1 cell array with matrix Qk per cell
% H: M x N cell array with channel Hk per cell
% order: 1 x K vector with BC encoding order
% Output
% S: K x 1 cell array of covariance matrix Sk per cell

% parameters
[M,N] = size(H{1});
K = length(H);

% initialization
S = cell(K,1);

I_M = eye(M);
I_N = eye(N);

pre = [];
S_sum = 0;

for k = K:-1:1
    current = find(order==k);
    x_index = find(order<k);
    Xk = I_N;
 
    if ~isempty(x_index)
        for i = x_index
            Xk = Xk + H{i}'*Q{i}*H{i};
        end
    end
    
    if ~isempty(pre)
        S_sum = S_sum + S{pre};
        Fk = I_M + H{current}*S_sum*H{current}';
        Hkeff = (Fk^(-1/2))'*H{current}*(Xk^(-1/2))';
        Qkeff = Fk^(1/2)*Q{current}*(Fk^(1/2))';
        S{current} = (Xk^(-1/2))'*ptpTransform(Qkeff,Hkeff)*Xk^(-1/2);
    else
        Fk = I_M + H{current}*S_sum*H{current}';
        Hkeff = (Fk^(-1/2))'*H{current}*(Xk^(-1/2))';
        Qkeff = Fk^(1/2)*Q{current}*(Fk^(1/2))';
        S{current} = (Xk^(-1/2))'*ptpTransform(Qkeff,Hkeff)*Xk^(-1/2);
    end
    
    pre = current;
    
end

