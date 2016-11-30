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
     (-k/m_t-k/m_p) -k/m_p -r*k/m_t (-c/m_t-c/m_p) -c/m_p -c*r/m_t 0;
     -k/m_p (-k/m_t-k/m_p) r*k/m_t -c/m_p (-c/m_t-c/m_p) c*r/m_t 0;
     -k*r/I k*r/I -2*(r^2)*k/I -c*r/I c*r/I -2*(r^2)*c/I 0;
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
rank(P);

Q = vertcat(C, C*A, C*(A^2), C(1,:)*(A^3));
rank(Q);

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

% calculate K
syms s

pole_1 = -.5;
pole_2 = -1;
pole_3 = -1.5;
pole_4 = -2;
pole_5 = -1;
pole_6 = -1.5;
pole_7 = -2;

expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4)*(s - pole_5)*(s - pole_6)*(s - pole_7))

a6 = 19/2;
a5 = 151/4;
a4 = 649/8;
a3 = 203/2;
a2 = 589/8;
a1 = 57/2;
a0 = 9/2;

% A_eq = [0 1 0 0 0 0 0;
%         0 0 1 0 0 0 0;
%         0 0 0 1 0 0 0;
%         0 0 0 0 1 0 0;
%         0 0 0 0 0 1 0;
%         0 0 0 0 0 0 1;
%         -a0 -a1 -a2 -a3 -a4 -a5 -a6;];
    
expand ((s - pole_1)*(s - pole_2)*(s - pole_3)*(s - pole_4))
expand ((s - pole_5)*(s - pole_6)*(s - pole_7))

a3 = 5;
a2 = 35/4;
a1 = 25/4;
a0 = 3/2;

a6 = 9/2;
a5 = 13/2;
a4 = 3;

A_eq = [0 1 0 0 0 0 0;
        0 0 1 0 0 0 0;
        0 0 0 1 0 0 0;
        -a0 -a1 -a2 -a3 0 0 0;
        0 0 0 0 0 1 0;
        0 0 0 0 0 0 1;
        0 0 0 0 -a4 -a5 -a6];

BK_co = A_co - A_eq

k11 = BK_co(4, 1) + BK_co(end, 1);
k12 = BK_co(4, 2) + BK_co(end, 2);
k13 = BK_co(4, 3) + BK_co(end, 3);
k14 = BK_co(4, 4) + BK_co(end, 4);
k15 = BK_co(4, 5) + BK_co(end, 5);
k16 = BK_co(4, 6) + BK_co(end, 6);
k17 = BK_co(4, 7) + BK_co(end, 7);

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


p = [pole_1 pole_2 pole_3 pole_4 pole_5 pole_6 pole_7];
K_matlab = place(A, B, p)

K_matlab_co = place(A_co, B_co, p)*T

BK = B*K
BK_matlab = B*K_matlab


% % Observer
% 
% % transform to observable canonical form
% U = [C(1,:); C(1,:)*A; C(1,:)*(A^2); C(1,:)*(A^3); C(2,:); C(2,:)*A; C(2,:)*(A^2)];
% U_inv = inv(U);
% 
% u1 = U_inv(:, 4);
% u2 = U_inv(:, end);
% 
% T_ob = horzcat(u1, A*u1, (A^2)*u1, (A^3)*u1, u2, A*u2, (A^2)*u2);
% T_ob_inv = inv(T_ob);
% 
% A_ob = T_ob_inv*A*T_ob
% transpose(A_co);
% 
% B_ob = T_ob_inv*B
% transpose(C_co);
% 
% C_ob = C*T_ob
% transpose(B_co)
% 
% % calculate L
% %A_eq_ob = A_ob - L*C_ob
% 
% pole_1_ob = 3*pole_1;
% pole_2_ob = 3*pole_2;
% pole_3_ob = 3*pole_3;
% pole_4_ob = 3*pole_4;
% pole_5_ob = 3*pole_5;
% pole_6_ob = 3*pole_6;
% pole_7_ob = 3*pole_7;
% 
% expand ((s - pole_1_ob)*(s - pole_2_ob)*(s - pole_3_ob)*(s - pole_4_ob)*(s - pole_5_ob)*(s - pole_6_ob)*(s - pole_7_ob))
% 
% a6_ob = 57/2;
% a5_ob = 1359/4;
% a4_ob = 17523/8;
% a3_ob = 16443/2;
% a2_ob = 143127/8;
% a1_ob = 41553/2;
% a0_ob = 19683/2;
% 
% % A_eq_ob = [0 0 0 0 0 0 -a0_ob;
% %             1 0 0 0 0 0 -a1_ob;
% %             0 1 0 0 0 0 -a2_ob;
% %             0 0 1 0 0 0 -a3_ob;
% %             0 0 0 1 0 0 -a4_ob;
% %             0 0 0 0 1 0 -a5_ob;
% %             0 0 0 0 0 1 -a6_ob];
%     
% expand ((s - pole_1_ob)*(s - pole_2_ob)*(s - pole_3_ob)*(s - pole_4_ob))
% expand ((s - pole_5_ob)*(s - pole_6_ob)*(s - pole_7_ob))
% 
% a3_ob = 15;
% a2_ob = 315/4;
% a1_ob = 675/4;
% a0_ob = 243/2;
% 
% a6_ob = 27/2;
% a5_ob = 117/2;
% a4_ob = 81;
% 
% A_eq_ob = [0 0 0 -a0_ob 0 0 0;
%             1 0 0 -a1_ob 0 0 0;
%             0 1 0 -a2_ob 0 0 0;
%             0 0 1 -a3_ob 0 0 0;
%             0 0 0 0 0 0 -a4_ob;
%             0 0 0 0 1 0 -a5_ob;
%             0 0 0 0 0 1 -a6_ob];
% 
% LC_ob = A_ob - A_eq_ob;
% L1 = LC_ob(:,4);
% L2 = LC_ob(:,end);
% 
% L_ob = horzcat(L1, L2);
% 
% L = T_ob*L_ob
% 
% p_ob = [pole_1_ob pole_2_ob pole_3_ob pole_4_ob pole_5_ob pole_6_ob pole_7_ob];
% L_matlab = place(A', C', p_ob).'
% 
% LC = L*C
