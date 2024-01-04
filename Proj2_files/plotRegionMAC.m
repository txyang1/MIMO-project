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
[~, I] = sort(R(1,:));
R_sort = R(:,I);
r_D = [0; R_sort(2,1)];
r_C = [R_sort(1,end); 0];
R_all = [r_D, R_sort, r_C];

plot(R_all(1,:),R_all(2,:),'Marker','square');
