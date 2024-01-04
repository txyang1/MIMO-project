clear all
close all
clc

% Load channel matrices from exampleMAC.mat
load('exampleMAC.mat');


% Transmit powers in dB
P_dB = [-10, 0, 10];

% Convert dB to linear scale
P = P.*10.^(P_dB / 10);

% Number of Pareto boundary sample points
S = 20;
exampleMACregion = figure(1);
hold on

  
%%
% Loop over different transmit powers
for i = 1:length(P)
    % Calculate Pareto boundary points
    R = ParetoBound(H, P(:,i), S);
    
    % Plot the capacity region
    exampleMACregion = plotRegionMAC(R,exampleMACregion);
    
    % Plot the Pareto boundary points
    hold on;
    
   
end

% Set axis labels and legend
xlabel('Rate User 1 (bps/Hz)');
ylabel('Rate User 2 (bps/Hz)');
legend({'P1 = -10 dB', 'P1 = 0 dB', 'P1 = 10 dB', 'Pareto Boundary'}, 'Location', 'Best');
% Save the figure
saveas(gcf, 'exampleMACregion.fig');
