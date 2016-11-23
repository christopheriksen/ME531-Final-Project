k = 500;  % spring coefficient
m = 70.0;  % mass of pod
c = 2*sqrt(k/m);  % damping coefficient
d = .5;  % diameter of pod
I = m*((d/2)^2)/2;  % inertia of pod

A = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 1 0 0;
    k/m k/m -2*k/m -2*c/m 0 0;
    0 0 0 0 0 1;
    -d*k/(2*I) d*k/(2*I) 0 0 0 0];

B = [1 0;
    0 1;
    0 0;
    c/m c/m;
    0 0
    -d*c/(2*I) d*c/(2*I)];

% determine if system is controllable
P = horzcat(B, A*B, (A^2)*B);
rank(P)

% transform to canonical form
M = horzcat(B(:, 1), A*(B(:, 1)), (A^2)*(B(:, 1)), B(:, 2), A*(B(:, 2)), (A^2)*(B(:, 2)));
M_inv = inv(M);

m1 = M_inv(3,:);
m2 = M_inv(end,:);

T = [m1; m1*A; m1*(A^2); m2; m2*A; m2*(A^2)];
T_inv = inv(T);

A_new = T*A*T_inv
B_new = T*B
%C_new = C*T_inv

syms s
pole_1 = -100;
pole_2 = -100;
pole_3 = -100;
pole_4 = -100;
pole_5 = -100;
pole_6 = -100;

expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6))

a5 = 8;
a4 = 28;
a3 = 56;
a2 = 68;
a1 = 48;
a0 = 16;

A_eq = [0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        -a0 -a1 -a2 -a3 -a4 -a5;];

% A_eq = A_new + B_new * K_new
%syms k11 k12 k13 k14 k15 k16 k21 k22 k23 k24 k25 k26
%K_new = [k11 k12 k13 k14 k15 k16;
%         k21 k22 k23 k24 k25 k26];
     
%BK = [0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      k11 k12 k13 k14 k15 k16;
%      0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      k21 k22 k23 k24 k25 k26]

BK = A_eq - A_new

k11 = BK(3, 1);
k12 = BK(3, 2);
k13 = BK(3, 3);
k14 = BK(3, 4);
k15 = BK(3, 5);
k16 = BK(3, 6);

k21 = BK(6, 1);
k22 = BK(6, 2);
k23 = BK(6, 3);
k24 = BK(6, 4);
k25 = BK(6, 5);
k26 = BK(6, 6);

K_new = [k11 k12 k13 k14 k15 k16;
         k21 k22 k23 k24 k25 k26];
     
K = K_new*T

B*K
