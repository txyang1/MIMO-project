function [Rwsum] = wsrQ(Q,H,w)
% Function [Rwsum] = wsrQ(Q,H,w)
%
% The function ratesMAC computes the rate region coordinates of
% the MIMO MAC for given transmit covariance matrices.
%
% Inputs
% Q: M x N x 2 array of the given tranmit covariance matrices
% H: M x N x 2 array of the given channel realizations
% w: 2 x 1 vector of weights
% Outputs
% Rwsum: achieveable weighted sum rate


M = size(H,1);

%TODO
Q1 = Q(:,:,1);
Q2 = Q(:,:,2);
H1 = H(:,:,1);
H2 = H(:,:,2);
if w(1) > w(2)
    Rwsum = real((w(1)-w(2))*log2(det(eye(M)+H1*Q1*H1'))+w(2)*log2(det(eye(M)+H1*Q1*H1'+H2*Q2*H2')));
elseif w(2) > w(1)
    Rwsum = real((w(2)-w(1))*log2(det(eye(M)+H2*Q2*H2'))+w(1)*log2(det(eye(M)+H1*Q1*H1'+H2*Q2*H2')));
end

%Team members: Tingxin Yang, Tian Yu
