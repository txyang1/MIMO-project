function [ Q,Cwsr ] = maxWSRmac(H,P,w)
%
% The function calculates the WSR capacity Cwsr and
% the optimal transmit covariance matrices Q of the
% two-user MAC with the convex optimization tool YALMIP
%
% Inputs
% H: M x N x 2 array of given channels to the users
% P: 2 x 1 column vector of available transmit power
% w: 2 x 1 column vector of weights w_1,...,w_K
% Outputs
% Cwsr: maximum WSR after convergence (WSR capacity)
% Q: N x N x 2 array of optimal covariance matrices
%
% Hint: Distinguish the cases w(2)>w(1) and w(1)>=w(2)

H1 = H(:,:,1);
H2 = H(:,:,2);
[M,N] = size(H1);

if w(2)>w(1)
    k = 2;
    j = 1;
else
    k = 1;
    j = 2;
end

Hk = H(:,:,k);
Hj = H(:,:,j);

options = sdpsettings('solver','sdpt3','verbose',0);

Qk = sdpvar(N,N,'hermitian','complex');
Qj = sdpvar(N,N,'hermitian','complex');

Constraints = [Qk>=0, Qj>=0, trace(Qk)<=P(k), trace(Qj)<=P(j)];
Objective = - (w(k)-w(j))*logdet(eye(M)+Hk*Qk*Hk') - w(j)*logdet(eye(M)+Hk*Qk*Hk'+Hj*Qj*Hj');

sol = optimize(Constraints, Objective, options);

Cwsr = real((w(k)-w(j))*log2(det(eye(M)+Hk*value(Qk)*Hk')) + w(j)*log2(det(eye(M)+Hk*value(Qk)*Hk'+Hj*value(Qj)*Hj')));

Q = zeros(N,N,2);
Q(:,:,k) = value(Qk);
Q(:,:,j) = value(Qj);

end

