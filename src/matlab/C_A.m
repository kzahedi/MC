function [ r ] = C_A( data, s_cardinality )
%C_A Summary of this function goes here
%   CA = ? p(s,a) p(s'|do(a)) log (p(s'|do(a)) / p(s'|do(s)))

eps     = 0.00000001;

s_max   = max(data(:,1));
a_max   = max(data(:,2));

s2      = data(2:end,1);
s1      = data(1:end-1,1);
a1      = data(1:end-1,2);

ps2s1a1 = sample3d([s2 s1 a1]); % add possibility to use other estimators

ps1a1   = sum(ps2s1a1, 1);
ps1     = sum(sum(ps2s1a1, 1),3);

% p(s'|s,a) = p(s',s,a) / p(s,a)
ps2_c_s1a1 = zeros(size(ps2s1a1));
for s2_index = 1:s_max
    for s1_index = 1:s_max
        for a1_index = 1:a_max
            if abs(ps1a1(1, s1_index, a1_index)) > eps
                ps2_c_s1a1(s2_index, s1_index, a1_index) = ps2s1a1(s2_index, s1_index, a1_index) / ps1a1(1, s1_index, a1_index);
            end
        end
    end
end


% p(s'|s,a) = p(s',s,a) / p(s,a)
ps2_c_s1a1 = zeros(size(ps2s1a1));
for s2_index = 1:s_max
    for s1_index = 1:s_max
        for a1_index = 1:a_max
            if abs(ps1a1(1, s1_index, a1_index)) > eps
                ps2_c_s1a1(s2_index, s1_index, a1_index) = ps2s1a1(s2_index, s1_index, a1_index) / ps1a1(1, s1_index, a1_index);
            end
        end
    end
end

% p(s'| do(a)) = ? p(s'|s,a) * p(s)
ps2_do_a1 = zeros(s_max, a_max);
for s2_index = 1:s_max
    for a1_index = 1:a_max
        for s1_index = 1:s_max
            ps2_do_a1(s2_index, a1_index) = ps2_do_a1(s2_index, a1_index) + ps2_c_s1a1(s2_index, s1_index, a1_index) * ps1(1,s1_index, 1);
        end
    end
end


% p(a|s) = p(s,a) / p(s)
pa1_c_s1 = zeros(a_max, s_max);
for a1_index = 1:a_max
    for s1_index = 1:s_max
        if abs(ps1(s1_index)) > eps
            pa1_c_s1(a1_index, s1_index) = ps1a1(1, s1_index, a1_index) / ps1(1, s1_index, 1);
        end
    end
end

% p(s'|do(s)) = ? p(a|s) p(s'|do(a))
ps2_do_s1 = zeros(s_max, s_max);
for s2_index = 1:s_max
    for s1_index = 1:s_max
        for a1_index = 1:a_max
            ps2_do_s1(s2_index, s1_index) = ps2_do_s1(s2_index, s1_index) + pa1_c_s1(a1_index, s1_index) * ps2_do_a1(s2_index, a1_index);
        end
    end
end

% CA = ? p(s,a) p(s'|do(a)) log (p(s'|do(a)) / p(s'|do(s)))
  r = 0.0;
  for s2_index = 1:s_max
    for s1_index = 1:s_max
      for a1_index = 1:a_max
        if abs(ps1a1(1, s1_index, a1_index)) > eps && abs(ps2_do_a1(s2_index, a1_index)) > eps && abs(ps2_do_s1(s2_index, s1_index)) > eps
          r = r + ps1a1(1, s1_index, a1_index) * ps2_do_a1(s2_index, a1_index) * (log(ps2_do_a1(s2_index, a1_index)) - log(ps2_do_s1(s2_index, s1_index)));
        end
      end
    end
  end
  r = 1.0 - r / log(s_cardinality);
end