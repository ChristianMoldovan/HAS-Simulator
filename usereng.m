function [ E1, E2 ] = usereng( rate, ratio, type)
%USERENG Summary of this function goes here
%   Detailed explanation goes here


if strcmp(type,'live')
    p = [20.211019015740877 1.736645481220426 52.834164850097640 22.769252206257367];
    q = [27.638403289303250 0.195201559665924 10.835323616531998 0.034081783312447];
else if strcmp(type,'vod')
        p = [25.237147323028220 14.256127913965312 12.516107735886832 0.258755664524164];
        q = [4.271232368792374 0.543520857704962 25.900008590831515 0.033852678414994];
    else disp('type = vod or live');
    end
end

E1 = p(1) * exp(-p(2) * rate) + p(3) * exp(-p(4) * rate); % rate
E2 = q(1) * exp(-q(2) * ratio) + q(3) * exp(-q(4) * ratio); % ratio

% E = (E1 + E2) / 2;
end

