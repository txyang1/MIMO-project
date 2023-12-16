function [R] = ParetoBound(H,P,S)
% Function [R] = ParetoBound(H,P,S)
%
% The function calculates S coordinate points at the
% Pareto Boundary of the two-user MIMO MAC capacity region
% via varying the weights w_1 and w_2.
%
% Inputs
% H: M x N x 2 array of given channels to the users
% P: 2 x 1 column vector of available transmit power
% S: number of Pareto Boundary sample points S
% Outputs
% R: 2 x 2S matrix of rate region coordinates

M = size(H,1);
N = size(H,2);

R = zeros(2,2*S);

w1 = linspace(0,1,2*S-1);
w2 = 1 - w1;


for index = 2:(length(w1)-1)
    Q = maxWSRmac(H,P,[w1(index), w2(index)]');
    R_bound = ratesMAC(Q,H);
    if w1(index)<w2(index)
        R(:,index) = R_bound(:,1);
    else
        if w1(index) == w2(index)
            R(:,index:index+1) = R_bound;
        else
            R(:,index+1) = R_bound(:,2);
        end
    end
end

k = 2;
j = 1;
Xk = H(:,:,k)' * H(:,:,k);
[Qk, Ck] = ratemaxQk(Xk, P(k));
Xj = H(:,:,j)' * (eye(M) + H(:,:,k)*Qk*H(:,:,k)')^-1 * H(:,:,j);
[~, Cj] = ratemaxQk(Xj, P(j));
R(:,1) = [Cj, Ck]';

k = 1;
j = 2;
Xk = H(:,:,k)' * H(:,:,k);
[Qk, Ck] = ratemaxQk(Xk, P(k));
Xj = H(:,:,j)' * (eye(M) + H(:,:,k)*Qk*H(:,:,k)')^-1 * H(:,:,j);
[~, Cj] = ratemaxQk(Xj, P(j));
R(:,end) = [Ck, Cj]';

end
    