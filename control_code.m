clc
clear
k = 500;  % spring coefficient
m_t = 150.0; % mass of thruster
m_p = 70.0;  % mass of pod
c = 2*sqrt(k/m_p);  % damping coefficient
d = .5;  % diameter of pod
I = m_p*((d/2)^2)/2;  % inertia of pod

A = [0 1 0 0 0 0 0 0;
     -k/m_t -c/m_t 0 0 k/m_t c/m_t -d*k/(2*m_t) -d*c/(2*m_t);
    0 0 0 1 0 0 0 0;
    0 0 -k/m_t -c/m_t k/m_t c/m_t d*k/(2*m_t) d*c/(2*m_t);
    0 0 0 0 0 1 0 0;
    k/m_p c/m_p k/m_p c/m_p -2*k/m_p -2*c/m_p 0 0;
    0 0 0 0 0 0 0 1;
    -d*k/(2*I) -d*c/(2*I) d*k/(2*I) d*c/(2*I) 0 0 -(d^2)*k/2 -(d^2)*c/2];

B = [0 0;
    1/m_t -1/(2*m_t);
    0 0;
    1/m_t 1/(2*m_t);
    0 0;
    0 0;
    0 0;
    0 0];

C = [0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 1 0];

% determine if system is controllable and observable
P = horzcat(B, A*B, (A^2)*B, (A^3)*B);
rank(P)

Q = vertcat(C, C*A, C*(A^2), C*(A^3));
rank(Q)

% transform to canonical form
M = horzcat(B(:, 1), A*(B(:, 1)), (A^2)*(B(:, 1)), (A^3)*(B(:, 1)), B(:, 2), A*(B(:, 2)), (A^2)*(B(:, 2)), (A^3)*(B(:, 2)));
M_inv = inv(M);

m1 = M_inv(4,:);
m2 = M_inv(end,:);

T = [m1; m1*A; m1*(A^2); m1*(A^3); m2; m2*A; m2*(A^2); m2*(A^3)];
T_inv = inv(T);

A_new = T*A*T_inv
B_new = T*B
C_new = C*T_inv

A_obs = transpose(A_new)
B_obs = transpose(C_new)
C_obs = transpose(B_new)

syms s
pole_1 = -1;
pole_2 = -1;
pole_3 = -1;
pole_4 = -1;
pole_5 = -1;
pole_6 = -1;
pole_7 = -1;
pole_8 = -1;

expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6)*(s - pole_7)*(s - pole_8))

a7 = 8;
a6 = 28;
a5 = 56;
a4 = 70;
a3 = 56;
a2 = 28;
a1 = 8;
a0 = 1;

A_eq = [0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 1;
        -a0 -a1 -a2 -a3 -a4 -a5 -a6 -a7;];

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

k11 = BK(4, 1);
k12 = BK(4, 2);
k13 = BK(4, 3);
k14 = BK(4, 4);
k15 = BK(4, 5);
k16 = BK(4, 6);
k17 = BK(4, 7);
k18 = BK(4, 8);

k21 = BK(8, 1);
k22 = BK(8, 2);
k23 = BK(8, 3);
k24 = BK(8, 4);
k25 = BK(8, 5);
k26 = BK(8, 6);
k27 = BK(8, 7);
k28 = BK(8, 8);

K_new = [k11 k12 k13 k14 k15 k16 k17 k18;
         k21 k22 k23 k24 k25 k26 k27 k28];
     
K = K_new*T

BK = B*K

%A_eq_obs = A_obs - L*C_obs

pole_1 = -10;
pole_2 = -10;
pole_3 = -10;
pole_4 = -10;
pole_5 = -10;
pole_6 = -10;
pole_7 = -10;
pole_8 = -10;

expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6)*(s - pole_7)*(s - pole_8))

a7_obs = 80;
a6_obs = 2800;
a5_obs = 56000;
a4_obs = 700000;
a3_obs = 5600000;
a2_obs = 28000000;
a1_obs = 80000000;
a0_obs = 100000000;

A_eq_obs = [0 0 0 0 0 0 0 -a0_obs;
            1 0 0 0 0 0 0 -a1_obs;
            0 1 0 0 0 0 0 -a2_obs;
            0 0 1 0 0 0 0 -a3_obs;
            0 0 0 1 0 0 0 -a4_obs;
            0 0 0 0 1 0 0 -a5_obs;
            0 0 0 0 0 1 0 -a6_obs;
            0 0 0 0 0 0 1 -a7_obs;]

LC = A_obs - A_eq_obs
