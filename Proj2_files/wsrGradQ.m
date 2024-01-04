function [V] = wsrGradQ(Q,H,w)
% Function [V] = wsrGradQ(Q,H,w)
%
% The function ratesMAC computes the rate region coordinates of
% the MIMO MAC for given transmit covariance matrices.
%
% Inputs
% Q: N x N x 2 array of the given tranmit covariance matrices
% H: M x N x 2 array of the given channel realizations
% w: 2 x 1 vector of weights
% Outputs
% V: N x N x 2 array of derivatives V1,V2


M = size(H,1);
N = size(H,2);

% TODO
V = zeros(N,N,2);
H1 = H(:,:,1);
H2 = H(:,:,2);
Q1 = Q(:,:,1);
Q2 = Q(:,:,2);

if w(1)>w(2)
    V(:,:,1) = ((w(1)-w(2))*H1'*inv(eye(M)+H1*Q1*H1')*H1+w(2)*H1'*inv(eye(M)+H2*Q2*H2'+H1*Q1*H1')*H1)*log2(exp(1));
    V(:,:,2) = log2(exp(1))*w(2)*H2'*inv(eye(M)+H1*Q1*H1'+H2*Q2*H2')*H2;
else
    V(:,:,1) = ((w(2)-w(1))*H2'*inv(eye(M)+H2*Q2*H2')*H2+w(1)*H2'*inv(eye(M)+H1*Q1*H1'+H2*Q2*H2')*H2)*log2(exp(1));
    V(:,:,2) = log2(exp(1))*w(1)*H1'*inv(eye(M)+H2*Q2*H2'+H1*Q1*H1')*H1;
end



