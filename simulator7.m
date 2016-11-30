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

BK =    [0         0         0         0         0         0         0;
         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0;
   -8.7265   -4.6233   -1.2358    7.6640   -1.9236    1.9617    0.4200;
   -4.6233   -8.7265    1.2358   -1.9236    7.6640   -1.9617    0.4200;
         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0;];

x0 = [0.0;0.0;0.0;0.0;0.0;0.0;0.0];
t = [0.0:0.1:100];
target = [0;0;pi/6;0;0;0.0;5.0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts);
plot(t, y(:,7), t, y(:,3), t, y(:,4), t, y(:,5));
%, t, y(:,3)

function xdot = controller(t, x, target)
    xdot = [0;
            0;
            0;
            0;
            0;
            0;
            0];
    f = 0;
    [f1,f2] = F12(x, target);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = x(6);
    xdot(4) = -k*(1/m_t+1/m_p)*x(1) - k/m_p*x(2) - r*k/m_t*sin(x(3)) - c*(1/m_t+1/m_p)*x(4) - c/m_p*x(5) - r*c/m_t*x(6)*cos(x(6)) + f1;
    xdot(5) = -k/m_p*x(1) - k*(1/m_t+1/m_p)*x(2) + r*k/m_t*sin(x(3)) - c/m_p*x(4) - c*(1/m_t+1/m_p)*x(5) + r*c/m_t*x(6)*cos(x(6)) + f2;
    xdot(6) = -k*r/I*x(1)*cos(x(3)) + k*r/I*x(2)*cos(x(3)) - 2*k*r*r*sin(x(3))*cos(x(3)) - c*r/I*x(4)*cos(x(3)) + c*r/I*x(5)*cos(x(3)) - 2*c*r*r/I*cos(2*x(3));
    xdot(7) = k/m_p*x(1) + k/m_p*x(2) + c/m_p*x(4) + c/m_p*x(5);
    
    
    
    function Fd = F(t)
        if t > 10.2
            Fd = 0.0;
        elseif t >= 10.0
            Fd = 100.0;
        else
            Fd = 0.0;
        end
    end
    function [F_1,F_2] = F12(x,target)
        x_err = target - x;
        x_err(1) = x_err(1) + r*sin(x(3));
        x_err(2) = x_err(2) - r*sin(x(3));
        x_err(4) = x_err(4) + r*x(6)*cos(x(3));
        x_err(5) = x_err(5) - r*x(6)*cos(x(3));
        action = BK*x_err;
        F_1 = action(4);
        F_2 = action(5);
    end
end

end
