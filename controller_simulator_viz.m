function main(export, skip)
clc
clear
        % Default so that we dont skip or export if we dont want to
        if nargin < 2

            skip = [0 0];

            if nargin < 1

                export = [0 0];

            end

        end

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


% BK =    [0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0;
%    -8.7265   -4.6233   -1.2358    7.6640   -1.9236    1.9617    0.4200;
%    -4.6233   -8.7265    1.2358   -1.9236    7.6640   -1.9617    0.4200;
%          0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0];

t = [0.0:0.1:10];
theta = 0.3;
vel = 1.0;

target = [-r*sin(theta); r*sin(theta); theta; 0; 0; 0.0; vel];
x0 =     [0.0; 0.0; 0; 0; 0; 0.0; 0.0];
% x0 =     [-r*sin(theta);r*sin(theta);theta;0;0;0.0;vel];

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts)
% plot(t,y(:,7), t,y(:,3), t,y(:,1), t,y(:,2))

% Define frame_info for graph
frame_info = graphs_create_elements

% Designate animation function
frame_gen_function = @(frame_info,tau) animation_graphs(t,y,frame_info,tau); % frame function, defined below in file

% Declare timing
timing.duration = 10; % three second animation
timing.fps = 15;     % create frames for 15 fps animation
timing.pacing = @(y) softspace(0,1,y); % Use a soft start and end, using the included softstart function

% Declare a directory name in which to place files
destination = 'Controller Starting At Rest';

% Animate the movie
[frame_info, endframe]...
    = animation(frame_gen_function,frame_info,timing,destination,export(1),skip(1));

function h = graphs_create_elements
    h.f = figure(1);                            % Designate a figure for this animation
    clf(h.f)                                     % Clear this figure
    set(h.f,'color','w','InvertHardCopy','off')  % Set this figure to have a white background
                                             %  and to maintain color
                                             %  settings during printing

    h.ax = axes('Parent',h.f);                   % Create axes for the plot
    set(h.ax,'Xlim',[0 10],'Ylim',[-1 1]);   % Set the range of the plot
    set(h.ax,'Xtick',0:1:10,'YTick',-1:.1:1);   % Set the tick locations
    set(h.ax,'FontSize',20);                       % Set the axis font size
    xlabel(h.ax, 't (s)')							 % Label the axes
%     ylabel(h.ax, 'v (m/s)')
    set(h.ax,'Box','on')						 % put box all the way around the axes
    title(h.ax, 'Controller Starting At Rest')

    % Line element to be used to draw the path
    h.line1 = line(0,0,'Color',[243, 156, 18]/255,'linewidth',1,'Parent',h.ax);
    h.line2 = line(0,0,'Color',[0 0 255]/255,'linewidth',1,'Parent',h.ax);
    h.line3 = line(0,0,'Color',[255 0 0]/255,'linewidth',1,'Parent',h.ax);
    h.line4 = line(0,0,'Color',[39, 174, 96]/255,'linewidth',1,'Parent',h.ax);
end

function frame_info = animation_graphs(t,y,frame_info,tau)
        step = .1;
        step_end = t(end);
        cur_step = t(1:round(tau*length(t)),:);
        index = round(((cur_step)/step)+1);

        p_velocity = y(:,7);
        p_ang = y(:,3);
        t1_position = y(:,1);
        t2_position = y(:,2);
        set(frame_info.line2,'XData',t(index),'YData',p_velocity(index))
        set(frame_info.line1,'XData',t(index),'YData',p_ang(index));
        set(frame_info.line3,'XData',t(index),'YData',t1_position(index));
        set(frame_info.line4,'XData',t(index),'YData',t2_position(index));
	frame_info.printmethod = @(dest) print(frame_info.f,'-dpng','-r 150','-painters',dest);
end

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
%     xdot(4) = (-(k/m_t)*(x(1) + r*sin(x(3))) -(k/m_p)*(x(1) + x(2))) + (-(c/m_t)*(x(4) + r*x(6)*cos(x(3))) -(c/m_p)*(x(4) + x(5))); %+ f1;
%     xdot(5) = (-(k/m_t)*(x(2) - r*sin(x(3))) -(k/m_p)*(x(1) + x(2))) + (-(c/m_t)*(x(5) - r*x(6)*cos(x(3))) -(c/m_p)*(x(4) + x(5))); %+ f2;

    xdot(4) = (-k*(x(1)+x(2)) - c*(x(4)+x(5)))/m_p + (-k*x(1) - c*x(4) - r*(k*sin(x(3)) + c*x(6)*cos(x(3))))/m_t + f1;
    xdot(5) = (-k*(x(1)+x(2)) - c*(x(4)+x(5)))/m_p + (-k*x(2) - c*x(5) + r*(k*sin(x(3)) + c*x(6)*cos(x(3))))/m_t + f2;


%     xdot(6) = (k*(-x(1) + x(2) -2*r*sin(x(3))) + c*(-x(4) + x(5) -2*r*x(6)*cos(x(3))))*(r*cos(x(3))/I);

    xdot(6) = (k*(x(2) - x(1)) + c*(x(5) - x(4)) -2*r*(k*sin(x(3)) + c*x(6)*cos(x(3))))*r*cos(x(3))/ I;
    xdot(7) = (k/m_p)*(x(1) + x(2)) + (c/m_p)*(x(4) + x(5)); % + f/m_p;


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
