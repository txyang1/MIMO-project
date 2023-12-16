function [fig] = plotRegionMAC(R,fig)
%
% Plots the rate region of a two user MIMO MAC with R1
% as abscissa and R2 as ordinate.
%
% Inputs
% R: 2 x S matrix of rate region coordinates
% fig: figure handle for plotting the boundary into a given figure
% Outputs
% fig: figure handle to the figure of the plotted boundary

if nargin<2, fig = figure; end

% TODO
figure(fig);

[~, I] = sort(R(1,:));
R_sort = R(:,I);
rD = [0, R_sort(2,1)].';
rC = [R_sort(1,end), 0].';
R_extend = [rD, R_sort, rC];
plot(R_extend(1,:), R_extend(2,:), 'LineWidth', 2.0, 'Marker', 'o', 'MarkerFaceColor', 'black');

