function [ r ] = MC_W( data, w_cardinality )
%MC_W Summary of this function goes here
%   Detailed explanation goes here

w2      = data(2:end,1);   % w'
w1      = data(1:end-1,1); % w
a1      = data(1:end-1,2); % a

pw2a1w1 = sample3d([w2 a1 w1]); % add possibility to use other estimators

r = 1 - CMI(pw2a1w1) / log(w_cardinality);

end