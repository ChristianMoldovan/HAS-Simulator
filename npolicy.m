function [ L,N,QoE ] = npolicy( mu,lambda, d )
%NPOLICY Summary of this function goes here
%   Detailed explanation goes here

a = lambda ./ mu;
EB0 = 1 ./ (mu .* (1 - a));
EBN = d .* EB0;
N = 1 ./ EBN;

EIN = d ./ lambda;
L = EIN;
QoE = exp(-(0.15 .* L + 0.19) .* N);

end

