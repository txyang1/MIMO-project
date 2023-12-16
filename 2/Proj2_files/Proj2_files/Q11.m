clear all
close all
clc

load('exampleMAC.mat');

Ptx = [0, 10, 20];
exampleMACregion = figure(1);

hold on
for i = 1:length(Ptx)
    R = ParetoBound(H,P*10^(Ptx(i)/10),20);
    exampleMACregion = plotRegionMAC(R, exampleMACregion);
end

legend('-10dB', '0dB' , '10dB')
hold off