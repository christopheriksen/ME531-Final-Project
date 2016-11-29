function simulator
clc
clear

k = 500;  % spring coefficient
m_t = 150.0; % mass of thruster
m_p = 70.0;  % mass of pod
c = 2*sqrt(k/m_p);  % damping coefficient
d = 1.0;  % diameter of pod
dt = 0.25; % diameter of thruster
Ip = m_p*((d/2)^2)/2;  % inertia of pod
It = m_t*((dt/2)^2)/2; % inertia of thruster
r = sqrt((d/2)^2+(d -(2*m_t*(d)/(2*m_t+m_p)))^2); % distance of thruster from cog
I = 2*(It + m_t*r^2) + Ip + m_p*(2*m_t*(d)/(2*m_t+m_p))^2; % inertia of system

BK =    [0         0         0         0         0         0         0         0;
   -2.0178    3.9529   -6.9179    0.8587    9.0407   -4.3753   -0.2143   -0.8905;
         0         0         0         0         0         0         0         0;
   -6.9179    0.8587   -2.0178    3.9529    9.0407   -4.3753    0.2143    0.8905;
         0         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0         0;
   -5.1453   -3.2490    5.1453    3.2490   -0.0000    0.0000    0.4500    1.8702;];

x0 = [0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;];
t = [0.0:0.1:50];
target = [0;0;0;0;0;0.0;1.0;0.0;];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts);
plot(t, y(:,7))

%t, y(:,2),t, y(:,4),t, y(:,6),
function xdot = controller(t, x, target)
    xdot = [0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;];
    f = 0;
    f1 = F1(x, target);
    f2 = F2(x, target);
    xdot(1) = x(2);
    xdot(2) = -k/m_t*x(1) - c/m_t*x(2) + k/m_t*x(5) + c/m_t*x(6) + f1/m_t;
    xdot(3) = x(4);
    xdot(4) = -k/m_t*x(3) - c/m_t*x(4) + k/m_t*x(5) + c/m_t*x(6) + f2/m_t;
    xdot(5) = x(6);
    xdot(6) = k/m_p*x(1) + c/m_p*x(2) + k/m_p*x(3) + c/m_p*x(4) - 2*k/m_p*x(5) - 2*c/m_p*x(6) + f/m_p;
    xdot(7) = x(8);
    xdot(8) = d*(f2 - f1)/(2*I);
    function Fd = F(t)
        if t > 10.2
            Fd = 0.0;
        elseif t >= 10.0
            Fd = 5000.0;
        else
            Fd = 0.0;
        end
    end
    function F_1 = F1(x,target)
        x_err = target - x;
        x_err(1) = 0.0;
        x_err(2) = 0.0;
        x_err(3) = 0.0;
        x_err(4) = 0.0;
        x_err(5) = 0.0;
        x_err(6) = 0.0;
        x_err(8) = 0.0;
        action = -BK*x_err;
        F_1 = action(2);
    end
    function F_2 = F2(x, target)
        x_err = target - x;
        x_err(1) = 0.0;
        x_err(2) = 0.0;
        x_err(3) = 0.0;
        x_err(4) = 0.0;
        x_err(5) = 0.0;
        x_err(6) = 0.0;
        x_err(8) = 0.0;
        action = -BK*x_err;
        F_2 = action(4);
    end
end

end
