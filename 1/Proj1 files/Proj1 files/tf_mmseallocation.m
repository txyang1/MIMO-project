function [ R, MSE ] = tf_mmseallocation( H, Cn, Ptx )
%
% MSE optimal power allocation with transmit filter design.
%
% Inputs
% H: Transmission Channel
% Ptx: available sum transmit power Ptx
% Outputs
% R: Achievable rate
% MSE: Achievable MSE



% TODO
[B, N] = size(H);
I_N = eye(N);

T = (H'*H + trace(Cn)/Ptx * I_N)^(-1) * H';
g_2 = 1/Ptx * trace((H'*H + trace(Cn)/Ptx*I_N)^(-2) * (H'*H));

MSE = B - trace(H*T) - trace(T'*H') + trace(T'*(H'*H)*T + Cn*g_2);
R = log2(det(H'*Cn^(-1)*H*(T*T')/g_2 + I_N));

MSE = real(MSE);
R = real(R);
end

