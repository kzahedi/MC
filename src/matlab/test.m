clf()
filename = '';
% read data
data = csvread(filename);

% extract data for W and A
position     = data(:,2);
velocity     = data(:,3);
accelaration = data(:,4);
action       = data(:,10);
sensor       = data(:,8);

% extract min/max values (over all experiments, to ensure comparability)
p_min = min(position);
p_max = max(position);

v_min = min(velocity);
v_max = max(velocity);

al_min = min(accelaration);
al_max = max(accelaration);

a_min = min(action);
a_max = max(action);

s_min = min(sensor);
s_max = max(sensor);


% binning for W and A
w_bins = 10; % this will become w_bins^(nr of variables)
a_bins = 100;
s_bins = 100;

d_position     = binVector(position,     p_min,  p_max,  w_bins);
d_velocity     = binVector(velocity,     v_min,  v_max,  w_bins);
d_accelaration = binVector(accelaration, al_min, al_max, w_bins);
d_action       = binVector(action,       a_min,  a_max,  a_bins);
d_sensor       = binVector(sensor,       s_min,  s_max,  s_bins);

%
% Plot binning for each variable
%
subplot(3,4,1)
plot(position(1:1000))
title('position')
subplot(3,4,2)
plot(d_position(1:1000))
title('discretised position')

subplot(3,4,3)
plot(velocity(1:1000))
title('velocity')
subplot(3,4,4)
plot(d_velocity(1:1000))
title('discretised velocity')

subplot(3,4,5)
plot(accelaration(1:1000))
title('accelaration')
subplot(3,4,6)
plot(d_accelaration(1:1000))
title('discretised accelaration')

subplot(3,4,7)
plot(action(1:1000))
title('action')
subplot(3,4,8)
plot(d_action(1:1000))
title('discretised action')

subplot(3,4,9)
plot(sensor(1:1000))
title('sensor')
subplot(3,4,10)
plot(d_sensor(1:1000))
title('discretised sensor')

%
% Combine accelartion, position, and velocity to a unary random variable
%
w = combineAndRelabelBinnedMatrix([d_position, d_velocity, d_accelaration]);

subplot(3,4,[11,12])
plot(w(1:5000))
title('discretised world state')

wa = [w d_action];
sa = [d_sensor d_action];

disp(sprintf('MC_A:   %f', MC_A(wa, w_bins^3)));
disp(sprintf('MC_W:   %f', MC_W(wa, w_bins^3)));
disp(sprintf('ASOC_A: %f', ASOC_A(sa, w_bins^3)));
disp(sprintf('ASOC_W: %f', ASOC_W(sa, w_bins^3)));
disp(sprintf('C_A:    %f', C_A(sa, w_bins^3)));
disp(sprintf('C_W:    %f', C_W(sa, w_bins^3)));


% tic(); disp(sprintf('MC_A:   %f', MC_A(wa, w_bins^3)));   toc()
% tic(); disp(sprintf('MC_W:   %f', MC_W(wa, w_bins^3)));   toc()
% tic(); disp(sprintf('ASOC_A: %f', ASOC_A(sa, w_bins^3))); toc()
% tic(); disp(sprintf('ASOC_W: %f', ASOC_W(sa, w_bins^3))); toc()
% tic(); disp(sprintf('C_A:    %f', C_A(sa, w_bins^3)));    toc()
% tic(); disp(sprintf('C_W:    %f', C_W(sa, w_bins^3)));    toc()
