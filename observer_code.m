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

C = [0 0 0 0 1 0 1 0];

% determine if system is observable
Q = vertcat(C, C*A, C*(A^2), C*(A^3), C*(A^4), C*(A^5), C*(A^6), C*(A^7))

%Q(1,:)
%Q(2,:)
%Q(3,:)
%Q(4,:)
%Q(5,:)
%Q(6,:)
%Q(7,:)
%Q(8,:)

rank(Q)

Q_inv = inv(Q);
q = Q_inv(:,1);
M_inv = horzcat(q, A*q, (A^2)*q, (A^3)*q, (A^4)*q, (A^5)*q, (A^6)*q, (A^7)*q);
M = inv(M_inv);

A_obs = M*A*M_inv

B_obs = M*B

C_obs = C*M_inv