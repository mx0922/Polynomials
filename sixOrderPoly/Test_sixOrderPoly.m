clear; close all; clc

% ��ʱ��
T = 1.0;

% ��ʼʱ��(t=0): λ�á��ٶȡ����ٶ�
p0 = 0;     v0 = 0;     a0 = 0;
% �м�ʱ��(t=T/2): λ��
pm = 0.2;
% ����ʱ��(t=T): λ�á��ٶȡ����ٶ�
p1 = 0.5;   v1 = 0;     a1 = 0;

% ���ζ���ʽ��ϵ��
S = getPolyCoeff(T, p0, v0, a0, pm, p1, v1, a1);

% ϸ��ʱ�̸����λ�á��ٶȡ����ٶ�
dt = 0.05;
tSpan = 0:dt:T;

Nt = length(tSpan);
pos = zeros(Nt, 1);
vel = zeros(Nt, 1);
acc = zeros(Nt, 1);

for i = 1:Nt
    t = tSpan(i);
    % �õ���Ӧʱ�̵�λ�á��ٶȡ����ٶ�
    [pos(i), vel(i), acc(i)] = getSixOrderPoly(S, t);
end

% % ����matlab�ľ�����������ֱ�ӵõ����ι켣
% [pos, vel, acc] = getSixOrderPoly(S, tSpan');

figure(1001); clf;
subplot(3, 1, 1)
plot(tSpan, pos)
subplot(3, 1, 2)
plot(tSpan, vel)
subplot(3, 1, 3)
plot(tSpan, acc)