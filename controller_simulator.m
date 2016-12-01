function simulator7
clc
clear

k = 500;  % spring coefficient
m_t = 50.0; % mass of thruster
m_p = 70.0;  % mass of pod
c = 2*sqrt(k/m_p);  % damping coefficient
r = .25;
I = m_p*((r)^2)/2;  % inertia of pod
% State [z1 z2 t z1d z2d td z3d]
A = [0              0              0            1              0              0           0;
     0              0              0            0              1              0           0;
     0              0              0            0              0              1           0;
    -(k/m_p+k/m_t) -k/m_p         -r*k/m_t     -(c/m_p+c/m_t) -c/m_p         -r*c/m_t     0;
    -k/m_p         -(k/m_p+k/m_t)  r*k/m_t     -c/m_p         -(c/m_t+c/m_p)  r*c/m_t     0;
    -k*r/I          k*r/I         -2*(r^2)*k/I -c*r/I          c*r/I         -2*(r^2)*c/I 0;
     k/m_p          k/m_p          0            c/m_p          c/m_p          0           0];
 
BK =    [-0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000;
   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000;
    0.0000   -0.0000    0.0000   -0.0000   -0.0000   -0.0000    0.0000;
    2.5904  -15.9403    4.2207    9.5877   -3.8473    2.5553    0.4200;
    6.6936  -20.0435    6.6923    0.0000    5.7404   -1.3680    0.4200;
    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000;
   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000;];

t = [0.0:0.01:100];
theta = 0.5;
vel = 1.0;

target = [-r*sin(theta);r*sin(theta);theta;0;0;0.0;vel];
x0 =     [0;0;theta;0;0;0.0;0.0];

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts);
plot(t, y(:,7), t, y(:,3))

function xdot = controller(t, x, target)
    xdot = [0;
            0;
            0;
            0;
            0;
            0;
            0];
    f = F(t);
    [f1,f2] = F12(x, target);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = (-k*(x(1)+x(2)) - c*(x(4)+x(5)))/m_p + (-k*x(1) - c*x(4) - r*(k*sin(x(3)) + c*x(6)*cos(x(3))))/m_t + f1;
    xdot(5) = (-k*(x(1)+x(2)) - c*(x(4)+x(5)))/m_p + (-k*x(2) - c*x(5) + r*(k*sin(x(3)) + c*x(6)*cos(x(3))))/m_t + f2;
    xdot(6) = r*cos(x(3)) * (k*(x(2) - x(1)) + c*(x(5) - x(4)) -2*r*(k*sin(x(3)) + c*x(6)*cos(x(3)))) / I;
    xdot(7) = k*(x(1) + x(2))/m_p + c*(x(4) + x(5))/m_p;% + f/m_p;

  
    function Fd = F(t)
        if t>30.0
            Fd = sin(t*pi/5);
        elseif t > 15.2
            Fd = 0.0;
        elseif t >= 15.0
            Fd = 200.0;
        else
            Fd = 0.0;
        end
    end
    function [F_1,F_2] = F12(x,target)
        x_err = target - x;
        action = BK*x_err;
        F_1 = action(4);
        F_2 = action(5);
    end
end
end