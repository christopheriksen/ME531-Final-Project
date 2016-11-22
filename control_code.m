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

% A_eq = A_new + B_new * K_new
