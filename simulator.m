function simulator
clc
clear

k = 300;  % spring coefficient
m_t = 50.0; % mass of thruster
m_p = 70.0;  % mass of pod
c = 2*sqrt(k/m_p);  % damping coefficient
d = .5;  % diameter of pod
I = m_p*((d/2)^2)/2;  % inertia of pod

BK =    [0         0         0         0         0         0         0         0;
  -10.1931    4.7397    4.2857    0.0591    6.0824   -4.0721   -3.6416    1.0793;
         0         0         0         0         0         0         0         0;
    4.2857    0.0591  -10.1931    4.7397    6.0824   -4.0721    3.6416   -1.0793;
         0         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0         0;
         0         0         0         0         0         0         0         0;];

x0 = [0.0;1.0;0.0;1.0;0.0;1.0;0.0;0.0;];
t = [0.0:0.1:200];
target = [0;0;0;0;0;1.0;0.0;0;];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts);
plot(t, y(:,2),t, y(:,4),t, y(:,6),t, y(:,7))


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
    xdot(2) = -k/m_t*x(1) - c/m_t*x(2) + k/m_t*x(5) + c/m_t*x(6) - d*k/(2*m_t)*sin(x(7)) - d*c/(2*m_t)*x(8)*cos(x(7)) + f1/m_t;
    xdot(3) = x(4);
    xdot(4) = -k/m_t*x(3) - c/m_t*x(4) + k/m_t*x(5) + c/m_t*x(6) + d*k/(2*m_t)*sin(x(7)) + d*c/(2*m_t)*x(8)*cos(x(7)) + f2/m_t;
    xdot(5) = x(6);
    xdot(6) = k/m_p*x(1) + c/m_p*x(2) + k/m_p*x(3) + c/m_p*x(4) - 2*k/m_p*x(5) - 2*c/m_p*x(6) + f/m_p;
    xdot(7) = x(8);
    xdot(8) = -d*k/(2*I)*x(1) - d*c/(2*I)*x(2) + d*k/(2*I)*x(3) + d*c/(2*I)*x(4) - (d^2)*k/2*sin(x(7)) - (d^2)*c/2*x(8)*cos(x(7));
    function Fd = F(t)
        if t > 10.2
            Fd = 0.0;
        elseif t >= 10.0
            Fd = 100.0;
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
        %x_err(8) = 0.0;
        x_err(6)
        x_err(7)
        action = -BK*x_err
        F_1 = action(2);
    end
    function F_2 = F2(x, target)
        x_err = target - x;
        x_err(1) = 0.0;
        x_err(2) = 0.0;
        x_err(3) = 0.0;
        x_err(4) = 0.0;
        x_err(5) = 0.0;
        %x_err(8) = 0.0;
        action = -BK*x_err;
        F_2 = action(4);
    end
end

end
