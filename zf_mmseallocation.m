function [ psi ] = zf_mmseallocation( phi, Ptx )
%
% MSE optimal power allocation under a zero-forcing constraint.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Ptx: available sum transmit power Ptx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN

sum_denom = 0; %The term on the denominator
for j=1:numel(phi)
    sum_denom = sum_denom + 1/sqrt(phi(j));
end

psi = zeros(size(phi)); %ready to save data

for i=1:numel(phi)    
    psi(i) = Ptx*(1/sqrt(phi(i)))/sum_denom;
end

%Team members: Tian Yu, Tingxin Yang