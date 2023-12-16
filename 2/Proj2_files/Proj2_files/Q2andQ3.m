clear all
close all
clc

load('exampleMAC.mat');

Ptx = [-10, 0, 10];
exampleMACregionQ = figure(1);
hold on
for i = 1:length(Ptx)
    [R, ~] = ratesMAC(10^(Ptx(i)/10)*Q, H);
    exampleMACregionQ = plotRegionMAC(R, exampleMACregionQ);
    
end

axis([0 5 0 5]);
legend('-10dB', '0dB' , '10dB')
hold off