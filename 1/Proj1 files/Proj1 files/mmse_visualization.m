% Script           mmse_visualization
%***************************************************************
%
% Visualization of the achievable MMSE of the considered four
% power allocation strategies for a point-to-point MIMO system.
% Plots the achievable MMSE versus the transmit power Etx in dB.
%
%***************************************************************

% System Parameters
clear all
clc

N = 4;
Ptx_dB = -20:0.1:30;
no_Ptx = length(Ptx_dB);

% Transmit power calculation
Ptx = 10 .^(Ptx_dB/10);% TODO

% Channel and eigenmode coefficients
load('example_channels.mat','H','Cn');
phi = eig( H' * inv(Cn) * H);
phi = phi(:); % TODO

% Initialization of MMSE arrays
MMSE_waterfilling = zeros(1,no_Ptx);
MMSE_mmse         = zeros(1,no_Ptx);
MMSE_uniform      = zeros(1,no_Ptx);
MMSE_tf_mmse      = zeros(1,no_Ptx);

% Calculation of the achievable MMSEs
% TODO
for i = 1:no_Ptx
    [psi_water,~,K_water] = waterfilling(phi, Ptx(i));
    
    MMSE_waterfilling(i) = sum(1 ./ (1+psi_water.*phi));
    
    [psi_uniform, K_uniform] = uniform_rate(phi,Ptx(i));
    
    MMSE_uniform(i) = sum(1 ./ (1+psi_uniform.*phi));
    
    [psi_mmse, mu_mmse, K_mmse] = mmseallocation(phi,Ptx(i));
    
    MMSE_mmse(i) = sum(1 ./ (1+psi_mmse.*phi));
    
    [~, MMSE_tf_mmse(i)] = tf_mmseallocation(H, Cn, Ptx(i));
end

% Calculation of the the transmit powers (in dB) where the number of streams switches from K to K+1
% TODO

Ptx_switch_waterfilling = activeStreams_waterfilling(phi);
Ptx_switch_mmse = activeStreams_mmse(phi);

Ptx_switch_waterfilling_in_dB = 10 * log10(Ptx_switch_waterfilling);
Ptx_switch_mmse_in_dB = 10 * log10(Ptx_switch_mmse);

% Calculation of the corresponding MMSEs to the switching powers from K to K+1
% TODO
MMSE_switch_waterfilling = zeros(size(Ptx_switch_waterfilling));
MMSE_switch_mmse = zeros(size(Ptx_switch_mmse));

for i = 1:(N-1)
    [psi_water,~,K_water] = waterfilling(phi, Ptx_switch_waterfilling(i));
    
    MMSE_switch_waterfilling(i) = sum(1 ./ (1+psi_water.*phi));
    
    [psi_mmse,mu_mmse, K_mmse] = mmseallocation(phi, Ptx_switch_mmse(i));
    
    MMSE_switch_mmse(i) = sum(1 ./ (1+psi_mmse.*phi));
end

% Plotting the achievable MMSEs over Etx in dB in a semilogarithmic scale
mmse_figure = figure;
hold on;
semilogy(Ptx_dB,MMSE_waterfilling,'Color','k','LineStyle','-','Marker','none','LineWidth',2);
semilogy(Ptx_dB,MMSE_mmse,'Color','r','LineStyle','-','Marker','none','LineWidth',2);
semilogy(Ptx_dB,MMSE_uniform,'Color','m','LineStyle','-','Marker','none','LineWidth',2);
semilogy(Ptx_dB,MMSE_tf_mmse,'Color','b','LineStyle','-','Marker','none','LineWidth',2);
hold off;

% Change visualization of the plot
xlabel('Etx in [dB]');
ylabel('mean square error');
legend(gca,'waterfilling','MMSE allocation','uniform','TF-MMSE allocation','Location','NorthWest');
grid on;

% Marking the switching points in the plots in the form
% semilogy(VALUES,VALUES,'Color',COLOR,'LineStyle','none','Marker',MARKER,'LineWidth',2);
% Use the same values for COLOR as are used for plotting the lines.
% Use the following values for MARKER: 
%        'o' for waterfilling
%        's' for mmse
hold on;
% TODO
plot(Ptx_switch_waterfilling_in_dB, MMSE_switch_waterfilling,'Color','k','LineStyle','none','Marker','o','LineWidth',2);
plot(Ptx_switch_mmse_in_dB, MMSE_switch_mmse,'Color','r','LineStyle','none','Marker','s','LineWidth',2);
hold off;

