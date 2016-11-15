function [WR] = hanning(len)

WR = .5*(1 - cos(2*pi*(1:len)'/(len+1)));

