% The Script calculates the sum capacity of the MIMO BC and the maximum
% achievable rates with block diagonalization with equal power allocation
% and block diagonlization with optimzal power allocation and plots the
% result versus Ptx

load('exampleBlockDiagMIMOBC.mat')

% transmit powers
Ptx = -10:2:40;
P = 10.^(Ptx/10);
no_P = length(P);

% Rate calculations
Req = zeros(1,no_P);
Rop = zeros(1,no_P);

for no = 1:no_P
    Req(no) = BlockDiagBCEqualPower( H,C,P(no) );
    Rop(no) = BlockDiagBC( H,C,P(no) );
end

figure;
hold on;

plot(Ptx,Req, ...
     'b-','LineWidth',1.5, ...
     'DisplayName',['Equal Power Allocation with Block Diagonalization Precoding']);

plot(Ptx,Rop, ...
     'r-','LineWidth',1.5, ...
     'DisplayName',['Optimal Power Allocation with Block Diagonalization Precoding']);
 
hold off;
xlabel('Ptx in [dB]');
ylabel('R in [bits/channel use]');
legend('show','Location','NorthWest');
