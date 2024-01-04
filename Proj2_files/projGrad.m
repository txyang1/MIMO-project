function [xOpt,fOpt] = projGrad(fun,grad,proj,x,nIter)
% Function [xOpt,fOpt] = projGrad(fun,grad,proj,x,nIter,method)
% 
% Projected gradient method.
%
% Inputs
% fun: function handle for the objective function f(x)
% grad: function handle for the gradient of the objective function
% proj: function handle for the projection onto the constraint set
% x: initial value for x
% nIter: maximum number of iterations
% method: optional parameter for the step size control
% Outputs
% xOpt: resulting optimization variables
% fOpt: corresponding value of the objective function

if nargin<5, nIter = 1e3; end
epsilon = 1e-4;
it = 0;
dist = Inf;
options = optimset('fminsearch');

while(it<nIter&&dist>epsilon)
    it = it+1;
    xOld = x;
    v = grad(x);
    [s,fval,exitflag] = fminsearch(@(z) fun(proj(x + z*v)),1,options);
    x = proj(x + s*v);
    dist = sum(abs((x(:)-xOld(:))/s).^2);
end

xOpt = x;
fOpt = fval;

%Team members: Tingxin Yang, Tian Yu
