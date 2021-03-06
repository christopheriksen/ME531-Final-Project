
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
    m_t = 150.0; % mass of thruster
    m_p = 70.0;  % mass of pod
    c = 2*sqrt(k/m_p);  % damping coefficient
    d = .5;  % diameter of pod
    r = .25;
    I = m_p*((d/2)^2)/2;  % inertia of pod

    BK =    [0         0         0         0         0         0         0         0;
      -10.1931    4.7397    4.2857    0.0591    6.0824   -4.0721   -3.6416    1.0793;
             0         0         0         0         0         0         0         0;
        4.2857    0.0591  -10.1931    4.7397    6.0824   -4.0721    3.6416   -1.0793;
             0         0         0         0         0         0         0         0;
             0         0         0         0         0         0         0         0;
             0         0         0         0         0         0         0         0;
             0         0         0         0         0         0         0         0;];


frame_info = graphs_create_elements

    
    x0 = [0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0];
    t = [0.0:0.1:200];
    target = [0;0;0;0;0;1;0.0;0;];
    opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
    [t, y] = ode45(@(t,y) controller(t,y,target), t, x0, opts);
    t_sys = t;
    y_sys = y;
    % plot(t, y(:,2),t, y(:,4),t, y(:,6),t, y(:,7));

    
% % Designate animation function
% frame_gen_function = @(frame_info,tau) animation_graphs(t,y,frame_info,tau); % frame function, defined below in file
% 
% % Declare timing
% timing.duration = 10; % three second animation
% timing.fps = 15;     % create frames for 15 fps animation
% timing.pacing = @(y) softspace(0,1,y); % Use a soft start and end, using the included softstart function
% 
% % Declare a directory name in which to place files
% destination = 'system_graphs';
% 
% % Animate the movie
% [frame_info, endframe]...
%     = animation(frame_gen_function,frame_info,timing,destination,export(1),skip(1));    


frame_info = system_create_elements

% Designate animation function
frame_gen_function = @(frame_info,tau) animation_system(t_sys,y_sys,frame_info,tau); % frame function, defined below in file
 
% Declare timing
timing.duration = 5; % three second animation
timing.fps = 15;     % create frames for 15 fps animation
timing.pacing = @(y) softspace(0,1,y); % Use a soft start and end, using the included softstart function

% Declare a directory name in which to place files
destination = 'system_plots';

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
    set(h.ax,'Xlim',[0 200],'Ylim',[0 1]);   % Set the range of the plot
    set(h.ax,'Xtick',0:20:200,'YTick',0:.1:1);   % Set the tick locations
    set(h.ax,'FontSize',20);                       % Set the axis font size
    xlabel(h.ax, 't (s)')							 % Label the axes
    ylabel(h.ax, 'v (m/s)')
    set(h.ax,'Box','on')						 % put box all the way around the axes
    title(h.ax, 'System Graphs')
    
    % Line element to be used to draw the path
    h.line1 = line(0,0,'Color',[255 0 150]/255,'linewidth',1,'Parent',h.ax);
    h.line2 = line(0,0,'Color',[0 0 255]/255,'linewidth',1,'Parent',h.ax);
    h.line3 = line(0,0,'Color',[255 0 0]/255,'linewidth',1,'Parent',h.ax);
    h.line4 = line(0,0,'Color',[255 150 100]/255,'linewidth',1,'Parent',h.ax);
end   

function h = system_create_elements
    h.f = figure(2);                            % Designate a figure for this animation
    clf(h.f)                                     % Clear this figure
    set(h.f,'color','w','InvertHardCopy','off')  % Set this figure to have a white background
                                             %  and to maintain color
                                             %  settings during printing

    h.ax = axes('Parent',h.f);                   % Create axes for the plot
    set(h.ax,'Xlim',[-5 5],'Ylim',[-5 5]);   % Set the range of the plot
    set(h.ax,'Xtick',-5:1:5,'YTick',-5:1:5);   % Set the tick locations
    set(h.ax,'FontSize',20);                       % Set the axis font size
    xlabel(h.ax, 'x (m)')							 % Label the axes
    ylabel(h.ax, 'y (m)')
    set(h.ax,'Box','on')						 % put box all the way around the axes
    axis equal
    title(h.ax, 'System Sim')
    
    % Line element to be used to draw the path
    h.line1 = line(0,0,'Color',[235 14 30]/255,'linewidth',2,'Parent',h.ax);
    h.line2 = line(0,0,'Color',[0 0 255]/255,'linewidth',2,'Parent',h.ax);
    h.line3 = line(0,0,'Color',[0 255 0]/255,'linewidth',2,'Parent',h.ax);
    h.line4 = line(0,0,'Color',[255 0 0]/255,'linewidth',2,'Parent',h.ax);
end   


function frame_info = animation_graphs(t,y,frame_info,tau)
        step = .1;
        step_end = t(end);
        cur_step = t(1:round(tau*length(t)),:);
        index = round(((cur_step)/step)+1);
        
        t1_velocity = y(:,2);
        t2_velocity = y(:,4);
        t3_velocity = y(:,6);
        t3_ang = y(:,7);
        t1_velocity(index)
        set(frame_info.line1,'XData',t(index),'YData',t1_velocity(index));
        set(frame_info.line2,'XData',t(index),'YData',t2_velocity(index));
        set(frame_info.line3,'XData',t(index),'YData',t3_velocity(index));
        set(frame_info.line4,'XData',t(index),'YData',t3_ang(index));
	frame_info.printmethod = @(dest) print(frame_info.f,'-dpng','-r 150','-painters',dest);
end

function frame_info = animation_system(t,y,frame_info,tau) 
    step = .1;
    step_end = t(end);
    cur_step_system = t(1:round(tau*length(t)),:);
    index_system = round(((cur_step_system)/step)+1);
    dy_timestep = (y(:,4)*(t(end)-t(end-1)));
    total_distace = cumsum(dy_timestep);

    array_x = repmat(0,1,2001);
    dx_timestep = (array_x(1,:)*(t(end)-t(end-1)));

%     step = .1;
%     step_end = t(end);

      center_y_timestep = 2+total_distace(index_system);
      y_timestep_d = center_y_timestep(1:round(tau*length(center_y_timestep)),:);
      x_y = y_timestep_d
      
      in = length(index_system)
      center_x_timestep = 2+dx_timestep(index_system);
      c_t = length(center_x_timestep)
      x_timestep_d = center_x_timestep(1:round(tau*length(center_x_timestep)),:);

      
%       x_timestep_t = array_x(1:round(tau*length(array_x)),:)
%     index_system = round(((cur_step)/step)+1)
%     d_timestep = (y(:,4)*(t(end)-t(end-1)));
%     total_distace = cumsum(d_timestep);

    if ~isempty(t)


%         t1 = get_square(x_timestep_d(end),y_timestep_d(end),1,1,0);  
%         set(frame_info.line1,'XData',t1(1,:),'YData',t1(2,:));
%         t2 = get_square(-2,2,1,1,0);
%         set(frame_info.line2,'XData',t2(1,:),'YData',t2(2,:));
%         p = get_square(0,0,1,1,0);
%         set(frame_info.line3,'XData',p(1,:),'YData',p(2,:));
    end
    
%     
%     set(frame_info.line1,'XData',t,'YData',y(:,2));
%     set(frame_info.line2,'XData',t,'YData',y(:,4));
%     set(frame_info.line3,'XData',t,'YData',y(:,6));
%     set(frame_info.line4,'XData',t,'YData',y(:,7));
	frame_info.printmethod = @(dest) print(frame_info.f,'-dpng','-r 150','-painters',dest);
end

function [square] = get_square(centerX, centerY,width,height, orientation)
    %        x3,x4
    %        -----   
    %  y1,y2 - c - y3,y4
    %        -----  
    %        x1,x2
    a = centerX;
    b = centerY;
    w = width;
    h = height;
    theta = orientation;
    X = [-w/2 w/2 w/2 -w/2 -w/2];
    Y = [h/2 h/2 -h/2 -h/2 h/2];
    P = [X;Y];
    ct = cos(theta);
    st = sin(theta);
    R = [ct -st;st ct];
    P = R * P;
    X_values = P(1,:)+a;
    gen_pt_num = 50;
    top_array = (linspace(X_values(1),X_values(2),gen_pt_num));
    right_array = (linspace(X_values(2),X_values(3),gen_pt_num));
    bottom_array = (linspace(X_values(3),X_values(4),gen_pt_num));
    left_array = (linspace(X_values(4),X_values(5),gen_pt_num));
    tot_x_array = [top_array,right_array,bottom_array,left_array];
    Y_values = P(2,:)+b;
  
    top_array = (linspace(Y_values(1),Y_values(2),gen_pt_num));
    right_array = (linspace(Y_values(2),Y_values(3),gen_pt_num));
    bottom_array = (linspace(Y_values(3),Y_values(4),gen_pt_num));
    left_array = (linspace(Y_values(4),Y_values(5),gen_pt_num));
    tot_y_array = [top_array,right_array,bottom_array,left_array];
    tot_array = [tot_x_array;tot_y_array];
%     h=plot(P(1,:)+a,P(2,:)+b);
%     axis equal;
    square = [tot_x_array;tot_y_array];
  
end
 

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
        x_err(1) = x(1) - (x(5)-r*sin(x(7)));
        x_err(2) = x(2) - (x(6)-r*x(8)*cos(x(7)));
        x_err(3) = x(3) - (x(5)+r*sin(x(7)));
        x_err(4) = x(4) - (x(6)+r*x(8)*cos(x(7)));
        x_err(5) = ((x(5) - x(1)) + (x(5) - x(3)))/2;
        action = -BK*x_err;
        F_1 = action(2);
    end
    function F_2 = F2(x, target)
        x_err = target - x;
        x_err(1) = x(1) - x(5);
        x_err(2) = x(2) - x(6);
        x_err(3) = x(3) - x(5);
        x_err(4) = x(4) - x(6);
        x_err(5) = ((x(5) - x(1)) + (x(5) - x(3)))/2;
        action = -BK*x_err;
        F_2 = action(4);
    end
end

end
