% Script           rate_visualization
%***************************************************************
%
% Visualization of the achievable rates of the considered four
% power allocation strategies for a point-to-point MIMO system.
% Plots the achievable rates versus the transmit power Ptx in dB.
%
%***************************************************************

% System Parameters
N = 4;
Ptx_dB = -20:0.1:30;
no_Ptx = length(Ptx_dB);

% Transmit power calculation
Ptx = db2pow(Ptx_dB);

% Channel and eigenmode coefficients
load('example_channels.mat','H','Cn');
phi = eig(H'*inv(Cn)*H);
phi = sort(phi,'descend');

% Initialization of rate arrays
R_waterfilling = zeros(1,no_Ptx);
R_uniform      = zeros(1,no_Ptx);
R_zf_mmse      = zeros(1,no_Ptx);

% Calculation of the achievable rates
% Function handle for rate calculation
rate = @(psi) sum(log2(1+phi.*psi));

for i = 1:numel(Ptx) %calculate for each Ptx value
    [psi_waterfilling,~,~] = waterfilling(phi,Ptx(i));
    R_waterfilling(i) = rate(psi_waterfilling);

    [psi_uniform,~] = uniform_rate(phi,Ptx(i));
    R_uniform(i) = rate(psi_uniform);

    psi_zf_mmse = zf_mmseallocation(phi,Ptx(i));
    R_zf_mmse(i) = rate(psi_zf_mmse);
end

% Calculation of the the transmit powers (in dB) where the number of streams switches from K to K+1
streampower_waterfilling = activeStreams_waterfilling(phi);

streampower_uniform = [];
for i=1:numel(Ptx)-1
    [~,k1] = uniform_rate(phi,Ptx(i));
    [~,k2] = uniform_rate(phi,Ptx(i+1));
    if k2>k1
        streampower_uniform = [streampower_uniform,Ptx(i+1)];
    end
end

streampowerdB_waterfilling = pow2db(streampower_waterfilling);
streampowerdB_uniform = pow2db(streampower_uniform);

% Calculation of the corresponding rates to the switching powers from K to K+1
streamrate_waterfilling = zeros(1,length(streampower_waterfilling));
for i=1:length(streampower_waterfilling)
    streamrate_waterfilling(i) = rate(waterfilling(phi,streampower_waterfilling(i)));
end

streamrate_uniform = zeros(1,length(streampowerdB_uniform));
for i=1:length(streampowerdB_uniform)
    streamrate_uniform(i) = rate(uniform_rate(phi,streampower_uniform(i)));
end


% Plotting the achievable rates over Ptx in dB
rate_figure = figure;
hold on;
plot(Ptx_dB,R_waterfilling,'Color','k','LineStyle','-','Marker','none','LineWidth',2);
plot(Ptx_dB,R_uniform,'Color','m','LineStyle','-','Marker','none','LineWidth',2);
plot(Ptx_dB,R_zf_mmse,'Color','b','LineStyle','-','Marker','none','LineWidth',2);
hold off;

% Changing the visualization of the plot
xlabel('Ptx in [dB]');
ylabel('rate in [bits/channel usage]');
legend(gca,'waterfilling','uniform','ZF-MMSE allocation','Location','NorthWest');
grid on;

% Marking the switching points in the plots in the form
% plot(VALUES,VALUES,'Color',COLOR,'LineStyle','none','Marker',MARKER,'LineWidth',2);
% Use the same values for COLOR as are used for plotting the lines.
% Use the following values for MARKER: 
%        'o' for waterfilling
%        's' for uniform
hold on;
plot(streampowerdB_waterfilling,streamrate_waterfilling,'Color','k','LineStyle','none','Marker','o','LineWidth',2);
plot(streampowerdB_uniform,streamrate_uniform,'Color','m','LineStyle','none','Marker','s','LineWidth',2);
hold off;

%Team members: Tian Yu, Tingxin Yang