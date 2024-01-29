function [fig] = plotSumRateBC(H,Ptx,fig)
% Function [fig] = plotSumRateBC(H,Ptx,fig)
%
% The function plots the sum capacity Csum and the individual
% achievable rates of each user in the MIMO BC for the encoding orders
% (1,...,K) and (K,...,1).
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% Ptx: column vector of Ptx values in dB
% fig: figure handle (optional)
% Outputs
% fig: figure handle

% system dimensions
[M,N] = size(H{1});
K = length(H);

% transmit powers
Ptx = Ptx(:);
P = 10.^(Ptx/10);
no_P = length(P);

% Rate calculations
Csum = zeros(no_P,1);
R_ascend = zeros(no_P,K);
R_descend = zeros(no_P,K);

for no = 1:no_P
    [Q,Csum(no)] = DualMACSumRateMaximization(H,P(no));
    [S_ascend] = MACtoBCtransform(Q,H,(1:K));
    [S_descend] = MACtoBCtransform(Q,H,(K:-1:1));
    R_ascend(no,:) = MAC_BC_rates(H,Q,S_ascend,(1:K))';
    R_descend(no,:) = MAC_BC_rates(H,Q,S_descend,(K:-1:1))';
end

if nargin<3, fig = figure; end
figure(fig);
hold on;
plot(Ptx,Csum, ...
     'b-','LineWidth',1.5, ...
     'DisplayName',['Csum']);

for i = 1:K
    ascend = sprintf('Ascend order: user %d',i);
    plot(Ptx,R_ascend(:,i),'LineWidth',1.5,'DisplayName',ascend); 
    descend = sprintf('Descend order: user %d',i);
    plot(Ptx,R_descend(:,i),'LineWidth',1.5,'DisplayName',descend);
end
 
hold off;
xlabel('Ptx in [dB]');
ylabel('R in [bits/channel use]');
legend('show','Location','NorthWest');
