function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
% function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
%
% The function calculates the sum capacity and the corresponding
% transmit covariance matrices of the K user MIMO MAC with a joint sum
% transmit power constraint.
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% Ptx: joint transmit power
%
% Outputs
% Q: K x 1 cell array of users' optimal transmit covariance Qk per cell
% Csum: sum capacity of the MIMO MAC

K = size(H,1);
[M,N] = size(H);


options = sdpsettings('solver','sdpt3','verbose',0);
Q = cell(K,1);

for k = 1:K
    Q{k} = sdpvar(M,M,'hermitian','complex');
end

Constraints = [];
tracesum = 0;
Z = eye(N);
for k = 1:K
    Constraints = [Constraints, Q{k}>=0];
    tracesum = tracesum + trace(Q{k});
    Z = Z + H{k}'*Q{k}*H{k};
end

Constraints = [Constraints, tracesum<=Ptx];
Objective = - logdet(Z);
sol = optimize(Constraints, Objective, options);

% P = Q
% P{1}
% P{2}
% value(P{1})
% [~,D] = eig(value(P{1}))

Q = cellfun(@value,Q,'UniformOutput',false);

[S] = MACtoBCtransform(Q,H,(1:K));
[ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, (1:K) );
Csum = sum(R_MAC);


% X = eye(N);
% %X_det = det(X)
% for k = 1:K
%     Hk = H{k}
%     Qk = Q{k}
%     Tk = H{k}'*Q{k}*H{k}*100
%     [~,Dk]=eig(H{k}'*Q{k}*H{k}*100)
%     X = H{k}'*Q{k}*H{k}+X
%     [~,Dx] = eig(X)
% %    de = prod(diag(D))
% %    det(X)
% end
% 
% Csum = real(log2(det(X)));

