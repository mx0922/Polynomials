clear; close all; clc

syms a b c d e f 'real'
syms p0 v0 p1 pm v1 a0 a1 'real'
syms t T 'real'

% 6�ζ���ʽ��λ�á��ٶȡ����ٶ�
pos = a + b * t + c * t^2 + d * t^3 + e * t^4 + f * t^5;
vel = diff(pos, t);
acc = diff(pos, t, 2);

% t=0: λ�á��ٶȡ����ٶ�
t = 0;
eq_p0 = subs(pos) == p0;
eq_v0 = subs(vel) == v0;
eq_a0 = subs(acc) == a0;

% t=T: λ�á��ٶȡ����ٶ�
t = T;
eq_p1 = subs(pos) == p1;
eq_v1 = subs(vel) == v1;
eq_a1 = subs(acc) == a1;

% ������ʽϵ��
S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1, eq_a1], [a,b,c,d,e,f]);
      
a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;
f = S.f;

matlabFunction(...
    a, b, c, d, e, f, ...
    'file', 'autoGen_fiveOrderPolyCoeff.m', ...
    'vars', {T, p0, v0, a0, p1, v1, a1});