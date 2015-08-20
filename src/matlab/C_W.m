function [ r ] = C_W( data, s_cardinality )
%C_W Summary of this function goes here
%   Detailed explanation goes here

eps     = 0.00000001;

s_max   = max(data(:,1));
a_max   = max(data(:,2));
s2      = data(2:end,1);
s1      = data(1:end-1,1);
a1      = data(1:end-1,2);

ps2s1a1 = sample3d([s2 s1 a1]); % add possibility to use other estimators

ps1     = sum(sum(ps2s1a1, 1),3);
pa1     = sum(sum(ps2s1a1, 1),2);
ps2s1   = sum(ps2s1a1, 3);
ps1a1   = sum(ps2s1a1, 1);
ps2a1   = sum(ps2s1a1, 2);

% p(s'|s) = p(s',s) / p(s)
ps2_c_s1 = zeros(s_max, s_max);
for s2_index = 1:s_max
    for s1_index = 1:s_max
        if abs(ps1(s1_index)) > eps
            ps2_c_s1(s2_index, s1_index) = ps2s1(s2_index, s1_index) / ps1(s1_index);
        end
    end
end

% p(s'|a) = p(s',a) / p(a)
ps2_c_a1 = zeros(s_max, a_max);
for s2_index = 1:s_max
    for a1_index = 1:a_max
        if abs(pa1(a1_index)) > eps
            ps2_c_a1(s2_index, a1_index) = ps2a1(s2_index, 1, a1_index) / pa1(a1_index);
        end
    end
end

% p(a|s) = p(a,s) / p(s)
pa1_c_s1 = zeros(a_max, s_max);
for a1_index = 1:a_max
    for s1_index = 1:s_max
        if abs(ps1(1, s1_index, 1)) > eps
            pa1_c_s1(a1_index, s1_index) = ps1a1(1, s1_index, a1_index) / ps1(1, s1_index, 1);
        end
    end
end

% p^(s'|s) = ? p(s'|a) p(a|s)
phat_s2_c_s1 = zeros(s_max, s_max);
for s2_index = 1:s_max
    for s1_index = 1:s_max
        for a1_index = 1:a_max
            phat_s2_c_s1(s2_index, s1_index) = ps2_c_a1(s2_index, a1_index) * pa1_c_s1(a1_index, s1_index);
        end
    end
end

% C_W = ? p(s',s) log p(s'|s) / phat(s'|s)
r = 0;
for s1_index = 1:s_max
    for s2_index = 1:s_max
        if abs(ps2s1(s2_index, s1_index, 1)) > eps && abs(ps2_c_s1(s2_index, s1_index)) > eps && abs(phat_s2_c_s1(s2_index, s1_index)) > eps
            r = r + ps2s1(s2_index, s1_index, 1) * (log(ps2_c_s1(s2_index, s1_index)) - log(phat_s2_c_s1(s2_index, s1_index)));
        end
    end
end

r = r / log(s_cardinality);

end

