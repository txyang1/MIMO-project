% Script           rate_visualization
%***************************************************************
%
% Visualization of the achievable rates of the considered four
% power allocation strategies for a point-to-point MIMO system.
% Plots the achievable rates versus the transmit power Ptx in dB.
%
%***************************************************************

% System Parameters
clear all
clc

N = 4;
Ptx_dB = -20:0.1:30;
no_Ptx = length(Ptx_dB);

% Transmit power calculation
Ptx = 10 .^(Ptx_dB/10);

% Channel and eigenmode coefficients
load('example_channels.mat','H','Cn');
phi = eig( H' * inv(Cn) * H);
phi = phi(:);

% Initialization of rate arrays
R_waterfilling = zeros(1,no_Ptx);
R_uniform         = zeros(1,no_Ptx);
R_mmse         = zeros(1,no_Ptx);
R_tf_mmse      = zeros(1,no_Ptx);

% Calculation of the achievable rates
% TODO
for i = 1:no_Ptx
    [psi_water,~,~] = waterfilling(phi, Ptx(i));
    R_waterfilling(i) = sum(log2(1 + psi_water .* phi));
    
    [psi_uniform,~] = uniform_rate(phi,Ptx(i));
    R_uniform(i) = sum(log2(1 + psi_uniform .* phi));
    
    [psi_mmse,~,~] = mmseallocation(phi,Ptx(i));
    R_mmse(i) = sum(log2(1 + psi_mmse .* phi));
    
    [R_tf_mmse(i),~ ] = tf_mmseallocation(H, Cn, Ptx(i));
end
% Calculation of the the transmit powers (in dB) where the number of streams switches from K to K+1
% TODO

Ptx_switch_waterfilling = activeStreams_waterfilling(phi);
Ptx_switch_mmse = activeStreams_mmse(phi);

Ptx_switch_waterfilling_in_dB = 10 * log10(Ptx_switch_waterfilling);
Ptx_switch_mmse_in_dB = 10 * log10(Ptx_switch_mmse);

% Calculation of the corresponding rates to the switching powers from K to K+1
% TODO
R_switch_waterfilling = zeros(size(Ptx_switch_waterfilling));
R_switch_mmse = zeros(size(Ptx_switch_mmse));

for i = 1:(N-1)
    [psi_water,~,~] = waterfilling(phi, Ptx_switch_waterfilling(i));
    R_switch_waterfilling(i) = sum(log2(1 + psi_water .* phi));
    
    [psi_mmse,~,~] = mmseallocation(phi, Ptx_switch_mmse(i));
    R_switch_mmse(i) = sum(log2(1 + psi_mmse .* phi));
end
    
% Plotting the achievable rates over Ptx in dB
rate_figure = figure;
hold on;
plot(Ptx_dB,R_waterfilling,'Color','k','LineStyle','-','Marker','none','LineWidth',2);
plot(Ptx_dB,R_mmse,'Color','r','LineStyle','-','Marker','none','LineWidth',2);
plot(Ptx_dB,R_uniform,'Color','m','LineStyle','-','Marker','none','LineWidth',2);
plot(Ptx_dB,R_tf_mmse,'Color','b','LineStyle','-','Marker','none','LineWidth',2);
hold off;

% Changing the visualization of the plot
xlabel('Ptx in [dB]');
ylabel('rate in [bits/channel usage]');
legend(gca,'waterfilling','MMSE allocation','uniform','TF-MMSE allocation','Location','NorthWest');
grid on;

% Marking the switching points in the plots in the form
% plot(VALUES,VALUES,'Color',COLOR,'LineStyle','none','Marker',MARKER,'LineWidth',2);
% Use the same values for COLOR as are used for plotting the lines.
% Use the following values for MARKER: 
%        'o' for waterfilling
%        's' for MMSE
hold on;
plot(Ptx_switch_waterfilling_in_dB, R_switch_waterfilling,'Color','k','LineStyle','none','Marker','o','LineWidth',2);
plot(Ptx_switch_mmse_in_dB, R_switch_mmse,'Color','r','LineStyle','none','Marker','s','LineWidth',2);
% TODO
hold off;
