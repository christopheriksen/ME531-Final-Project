clc
clear
k = 500;  % spring coefficient
m_t = 50.0; % mass of thruster
m_p = 70.0;  % mass of pod
c = 2*sqrt(k/m_p);  % damping coefficient
r = .25;  % radius of pod
I = m_p*((r)^2)/2;  % inertia of pod

A = [0 0 0 1 0 0 0;
     0 0 0 0 1 0 0;
     0 0 0 0 0 1 0;
     -k/m_t 0 0 -c/m_t 0 0 0;
     0 -k/m_t 0 0 -c/m_t 0 0;
     -k*r/I k*r/I 0 -c*r/I c*r/I 0 0;
     k/m_p k/m_p 0 c/m_p c/m_p 0 0];

B = [0 0;
     0 0;
     0 0;
     1/m_t 0;
     0 1/m_t;
     0 0;
     0 0];
    

C = [0 0 1 0 0 0 0;
     0 0 0 0 0 0 1];

% determine if system is controllable and observable
P = horzcat(B, A*B, (A^2)*B, (A^3)*B(:,1));
rank(P)

Q = vertcat(C, C*A, C*(A^2), C(1,:)*(A^3));
rank(Q)

% transform to controllable canonical form
M = horzcat(B(:, 1), A*(B(:, 1)), (A^2)*(B(:, 1)), (A^3)*(B(:, 1)), B(:, 2), A*(B(:, 2)), (A^2)*(B(:, 2)));
M_inv = inv(M);

m1 = M_inv(4,:);
m2 = M_inv(end,:);

T = [m1; m1*A; m1*(A^2); m1*(A^3); m2; m2*A; m2*(A^2)];
T_inv = inv(T);

A_co = T*A*T_inv
B_co = T*B
C_co = C*T_inv

% transform to observable canonical form
U = [C(1,:); C(1,:)*A; C(1,:)*(A^2); C(1,:)*(A^3); C(2,:); C(2,:)*A; C(2,:)*(A^2)];
U_inv = inv(U);

u1 = U_inv(:, 4);
u2 = U_inv(:, end);

T_ob = horzcat(u1, A*u1, (A^2)*u1, (A^3)*u1, u2, A*u2, (A^2)*u2);
T_ob_inv = inv(T_ob);

A_obs = T_ob_inv*A*T_ob;
transpose(A_co);

B_obs = T_ob_inv*B;
transpose(C_co);

C_obs = C*T_ob;
transpose(B_co);


% calculate K
syms s

pole_1 = -.5;
pole_2 = -1;
pole_3 = -1.5;
pole_4 = -2;
pole_5 = -.5;
pole_6 = -1;
pole_7 = -1.5;

expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6)*(s - pole_7))

a6 = 8;
a5 = 53/2;
a4 = 47;
a3 = 769/16;
a2 = 113/4;
a1 = 141/16;
a0 = 9/8;

A_eq = [0 1 0 0 0 0 0;
        0 0 1 0 0 0 0;
        0 0 0 1 0 0 0;
        0 0 0 0 1 0 0;
        0 0 0 0 0 1 0;
        0 0 0 0 0 0 1;
        -a0 -a1 -a2 -a3 -a4 -a5 -a6;];
    
expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4))
expand ((s - pole_5)*(s - pole_6)*(s - pole_7))

a3 = 5;
a2 = 35/4;
a1 = 25/4;
a0 = 3/2;

a6 = 3;
a5 = 11/4;
a4 = 3/4;

A_eq = [0 1 0 0 0 0 0;
        0 0 1 0 0 0 0;
        0 0 0 1 0 0 0;
        -a0 -a1 -a2 -a3 0 0 0;
        0 0 0 0 0 1 0;
        0 0 0 0 0 0 1;
        0 0 0 0 -a4 -a5 -a6];

% A_eq = A_new + B_new * K_new
%syms k11 k12 k13 k14 k15 k16 k21 k22 k23 k24 k25 k26
%K = [k11 k12 k13 k14 k15 k16;
%         k21 k22 k23 k24 k25 k26];
     
%BK = [0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      k11 k12 k13 k14 k15 k16;
%      0 0 0 0 0 0;
%      0 0 0 0 0 0;
%      k21 k22 k23 k24 k25 k26]

BK_co = A_co - A_eq;

k11 = BK_co(4, 1);
k12 = BK_co(4, 2);
k13 = BK_co(4, 3);
k14 = BK_co(4, 4);
k15 = BK_co(4, 5);
k16 = BK_co(4, 6);
k17 = BK_co(4, 7);

k21 = BK_co(end, 1);
k22 = BK_co(end, 2);
k23 = BK_co(end, 3);
k24 = BK_co(end, 4);
k25 = BK_co(end, 5);
k26 = BK_co(end, 6);
k27 = BK_co(end, 7);

K_co = [k11 k12 k13 k14 k15 k16 k17;
         k21 k22 k23 k24 k25 k26 k27];
     
K = K_co*T


p = [-.5 -1 -1.5 -2 -1 -1.5 -2];
K_matlab = place(A, B, p)

BK = B*K
BK_matlab = B*K_matlab

% calculate L
%A_eq_obs = A_obs - L*C_obs

% pole_1 = -2;
% pole_2 = -2;
% pole_3 = -2;
% pole_4 = -2;
% pole_5 = -2;
% pole_6 = -2;
% pole_7 = -2;
% pole_8 = -2;
% 
% expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6)*(s - pole_7)*(s - pole_8))
% 
% a7_obs = 30;
% a6_obs = 765/2;
% a5_obs = 2700;
% a4_obs = 184113/16;
% a3_obs = 241785/8;
% a2_obs = 761805/16;
% a1_obs = 164025/4;
% a0_obs = 59049/4;
% 
% A_eq_obs = [0 0 0 0 0 0 0 -a0_obs;
%             1 0 0 0 0 0 0 -a1_obs;
%             0 1 0 0 0 0 0 -a2_obs;
%             0 0 1 0 0 0 0 -a3_obs;
%             0 0 0 1 0 0 0 -a4_obs;
%             0 0 0 0 1 0 0 -a5_obs;
%             0 0 0 0 0 1 0 -a6_obs;
%             0 0 0 0 0 0 1 -a7_obs];
%         
% expand ((s - pole_1)*(s - pole_3)*(s - pole_5)*(s - pole_7))
% 
% a3_obs = 15;
% a2_obs = 315/4;
% a1_obs = 675/4;
% a0_obs = 243/2;
% 
% A_eq_obs = [0 0 0 -a0_obs 0 0 0 0;
%             1 0 0 -a1_obs 0 0 0 0;
%             0 1 0 -a2_obs 0 0 0 0;
%             0 0 1 -a3_obs 0 0 0 0;
%             0 0 0 0 0 0 0 -a0_obs;
%             0 0 0 0 1 0 0 -a1_obs;
%             0 0 0 0 0 1 0 -a2_obs;
%             0 0 0 0 0 0 1 -a3_obs];
% 
% LC_obs = A_obs - A_eq_obs;
% L1 = LC_obs(:,4);
% L2 = LC_obs(:,8);
% 
% L_obs = horzcat(L1, L2);
% 
% L = T_ob*L_obs
% 
% pl = [-1.5 -1.5 -3 -3 -4.5 -4.5 -6 -6];
% %L_matlab = place(A', C', pl).'
% 
% LC = L*C
% 
% 
% % controller-observer
% %A_eq = A + BK - LC
% %LC_obs = A_obs + T_ob_inv*B*K*T_ob - A_eq_obs
% %L1 = LC_obs(:,4);
% %L2 = LC_obs(:,8);
% 
% %L_obs = horzcat(L1, L2);
% 
% %L = T_ob_inv*L_obs;
% %L1 = L(:,1)
% %L2 = L(:,2)
% 
% %LC = L*C
