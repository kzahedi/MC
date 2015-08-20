function [ r ] = MC_A( data, w_cardinality )
%MC_A Summary of this function goes here
%   Detailed explanation goes here

w2      = data(2:end,1);  % w;
w1      = data(1:end-1,1); % w
a1      = data(1:end-1,2); % a

pw2w1a1 = sample3d([w2 w1 a1]); % add possibility to use other estimators

r = CMI(pw2w1a1) / log(w_cardinality);

end

