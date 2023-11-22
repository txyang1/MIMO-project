function [ psi ] = zf_mmseallocation( phi, Ptx )
%
% MSE optimal power allocation under a zero-forcing constraint.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Ptx: available sum transmit power Ptx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
  N = length(phi);
  psi = zeros(size(phi));
  for i = 1:N
       for  j = 1:N
         psi(i) = Ptx/(sum(1/sqrt(phi(j)))*sqrt(phi(i)));
       end 
  end
       
%TODO

end

