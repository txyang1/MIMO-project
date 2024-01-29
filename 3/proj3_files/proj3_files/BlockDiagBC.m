function [ Rsum ] = BlockDiagBC( H,C,Ptx )
% function [ Rsum ] = BlockDiagBC( H,C,Ptx )
%
% The function computes the maximum achievable rate of a MIMO BC with block
% diagonalization.
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% C: K x 1 cell array with noise covariance C_nk per cell
% Ptx: joint available transmit power
% Outputs
% Rsum: achievable sum rate

[M,N] = size(H{1});
K = length(H);
phi = [];

for k = 1:K
    H_ = [];
    for j = 1:K
        if j~=k
            H_ = [H_;H{j}];
        else
            continue
        end
    end

    Pk = eye(N) - H_'*(H_*H_')^(-1)*H_;
    [Nk,Dk] = eig(Pk);
    dk = round(real(diag(Dk)));

    Vk = Nk(:,dk==1);
    Xk = Vk'*H{k}'*(C{k})^(-1)*H{k}*Vk;
    [~,Phik] = eig(Xk);
    phik = real(diag(Phik));
    phi = [phi;phik];
end

psi = waterfilling(phi,Ptx);

Rsum = real(sum(log2(1+phi.*psi),'all'));
% options = sdpsettings('solver','sdpt3','verbose',0);
% 
% Q = cell(K,1);
% for k = 1:K
%     Q{k} = sdpvar(M,M,'hermitian','complex');
% end
% 
% Constraints = [];
% tracesum = 0;
% 
% for k = 1:K
%     Constraints = [Constraints, Q{k}>=0];
%     tracesum = tracesum + trace(Q{k});
% end
% Constraints = [Constraints, tracesum<=Ptx];
% 
% Objective = 0;
% for k = 1:K
%     Objective = Objective - logdet(eye(M)+X{k}*Q{k});
% end
% sol = optimize(Constraints, Objective, options);
% 
% Q = cellfun(@value,Q,'UniformOutput',false);
% 
% Rsum = 0;
% for k = 1:K
%     Rsum = Rsum + log2(det(eye(M) + X{k}*Q{k}));
% end










