clear; close all; clc

% total time
T = 1.0;

% start time(t=0): pos, vel, acc
p0 = 0;     v0 = 0;     a0 = 0;
% middle time(t=T/2):pos
pm = 0.2;
% end time(t=T):pos, vel, acc
p1 = 0.5;   v1 = 0;     a1 = 0;

% polynomials coefficients
S = getPolyCoeff(T, p0, v0, a0, pm, p1, v1, a1);

% finer time node
dt = 0.05;
tSpan = 0:dt:T;

Nt = length(tSpan);
pos = zeros(Nt, 1);
vel = zeros(Nt, 1);
acc = zeros(Nt, 1);

for i = 1:Nt
    t = tSpan(i);
    % finer time node pos, vel, acc
    [pos(i), vel(i), acc(i)] = getSixOrderPoly(S, t);
end

% % good vectorization
% [pos, vel, acc] = getSixOrderPoly(S, tSpan');

figure(1001); clf;
subplot(3, 1, 1)
plot(tSpan, pos)
subplot(3, 1, 2)
plot(tSpan, vel)
subplot(3, 1, 3)
plot(tSpan, acc)