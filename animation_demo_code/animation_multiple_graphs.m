function animation_demonstration(export,skip)
% Produces two simple movies, to demonstrate usage of the ANIMATION
% function for exporting movie frames
%
%
% ========================
% Copyright Ross L. Hatton, 2011


	% Handle default arguments -- don't export frames, and don't skip
	% unless told otherwise
	if nargin < 2
		
		skip = [0 0];
		
		if nargin < 1
			
			export = [0 0];
			
		end
		
	end

	%%%%%%%%
	% Create animation elements, and store them in the frame_info structure
	frame_info = animation_demo_create_elements; %setup function, defined below in file

	
	%%%%%%%%%%%%%%%%%
	% First movie: Draw the cosine function
	
	% Designate animation function
	frame_gen_function...
		= @animation_demo_draw_cosine; % frame function, defined below in file

	% Declare timing
	timing.duration = 5; % three second animation
	timing.fps = 15;     % create frames for 15 fps animation
	timing.pacing = @(y) softspace(0,1,y); % Use a soft start and end, using the included softstart function

	% Declare a directory name in which to place files
	destination = 'demonstration_movie_1';
	
	% Animate the movie
	[frame_info, endframe]...
		= animation(frame_gen_function,frame_info,timing,destination,export(1),skip(1));
	
% 	%%%%%%%%%%%%%%%%
% 	% Second movie: Zoom in on first quadrant
% 	
% 	% Designate animation function. Anonymous function format allows for
% 	% specifing the starting limits based on the current state of the axes
% 	
% 	start_limits = axis(frame_info.ax); % The axis limits at the end of the previous movie
% 	
% 	frame_gen_function...
% 		= @(frame_info,tau) animation_demo_panzoom(start_limits, frame_info, tau);
% 
% 	% Declare timing
% 	timing.duration = 3; % three second animation
% 	timing.fps = 15;     % create frames for 30 fps animation
% 	timing.pacing = @(y) softspace(0,1,y); % Use a soft start and end, using the included softstart function
% 
% 	% Declare a directory name in which to place files
% 	destination = 'demonstration_movie_2';
% 
% 	% Animte the movie
% 	[frame_info, endframe] = animation(frame_gen_function,frame_info,timing,destination,export(2),skip(2),endframe);
	
end


%%%%%%%%%%%%%%%%%%
% Create animation elements, and return a vector of their handles
function h = animation_demo_create_elements

	h.f = figure(17);                            % Designate a figure for this animation
	clf(h.f)                                     % Clear this figure
	set(h.f,'color','w','InvertHardCopy','off')  % Set this figure to have a white background
												 %  and to maintain color
												 %  settings during printing

	h.ax = axes('Parent',h.f);                   % Create axes for the plot
	set(h.ax,'Xlim',[0 10],'Ylim',[0 10]);   % Set the range of the plot
	set(h.ax,'Xtick',0:1:10,'YTick',0:1:10);   % Set the tick locations
	set(h.ax,'FontSize',20);                       % Set the axis font size
	xlabel(h.ax, 'x')							 % Label the axes
	ylabel(h.ax, 'y')
	set(h.ax,'Box','on')						 % put box all the way around the axes


	% Line element to be used to draw the cosine function
    %Line 1 Definition
	h.line = line(0,0,'Color',[235 14 30]/255,'linewidth',5,'Parent',h.ax);
    %Line 2 Definition
    h.line_2 = line(0,0,'Color',[0 0 235]/255,'linewidth',5,'Parent',h.ax);
    %Line 3 Definition
    h.line_3 = line(0,0,'Color',[0 0 235]/255,'linewidth',5,'Parent',h.ax);
    %Line 4 Definition
    h.line_4 = line(0,0,'Color',[0 0 235]/255,'linewidth',5,'Parent',h.ax);
    %Line 5 Definition
    h.line_5 = line(0,0,'Color',[0 0 235]/255,'linewidth',5,'Parent',h.ax);

end


%%%%%%%%%%%%%%%%%%
% Frame content specification for dynamically drawing the cosine function,
% at time tau on a range of 0 to 1
function frame_info = animation_demo_draw_cosine(frame_info,tau)

	% Declare a baseline set of points at which to evaluate the cosine
	% function
	x_full = linspace(0,10,300)';
	
	% Truncate the points to the percentage of the way through drawing
	x = x_full(1:round(tau*length(x_full)),:);

	% Evaluate the cosine function at the x points;
	y = sin(x)+5;
	z = sin(x)
    %Combining array testing for ODE45 output
%     s = [y,z]
    
	% Set the line in the plot to the calculated x and y points
	%Plot 1
    set(frame_info.line,'XData',x,'YData',y);
	%Plot 2
%     set(frame_info.line_2,'XData',x,'YData',z);
    if ~isempty(x)
        square = get_square(x(end),y(end),0);
        set(frame_info.line_2,'XData',square(:,1),'YData',square(:,2))
	    set(frame_info.line_3,'XData',square(:,3),'YData',square(:,4))
        set(frame_info.line_4,'XData',square(:,5),'YData',square(:,6))
        set(frame_info.line_5,'XData',square(:,7),'YData',square(:,8))
    %Set the line in the plot to the two columns in s, this can be chaged
    %depending on input form.
% 	set(frame_info.line,'XData',x,'YData',s(:,1))
% 	set(frame_info.line_2,'XData',x,'YData',s(:,2))

	% Declare a print method (in this case, print 150dpi png files of
	% frame_info.f, using the Painters renderer)
	frame_info.printmethod = @(dest) print(frame_info.f,'-dpng','-r 150','-painters',dest);
	

end


% Function to produce an ellipse at the point where the vehicle is
% Size of the ellipse is predefined based upon how big the vehicle is

function [ellipse] = get_ellipse(centerX, centerY, orientation)
    
    xAxis = .75;
    yAxis = .11;
    
    theta = linspace(0, 2*pi, 150);
    
    xx = (xAxis) * sin(theta) + centerX;
    yy = (yAxis) * cos(theta) + centerY;
    
    xx2 = (xx-centerX)*cos(orientation) - (yy-centerY)*sin(orientation) + centerX;
    yy2 = (xx-centerX)*sin(orientation) + (yy-centerY)*cos(orientation) + centerY;
    
    ellipse = [xx2',yy2'];
    
end

function [square] = get_square(centerX, centerY, orientation)
    %        x3,x4
    %        -----   
    %  y1,y2 - c - y3,y4
    %        -----  
    %        x1,x2 
    xsize = .5;
    ysize = .5;
    
    line = linspace(0, .5, 10);
    
    x1 = (1)* (line) + centerX -.25
    x2 = (0)* (line) + centerY 
    y1 = (0)* (line) + centerX 
    y2 = (1)* (line) + centerY 
    x3 = (1)* (line) + centerX -.25
    x4 = (0)* (line) + centerY +.5 
    y3 = (0)* (line) + centerX +.5 -.75
    y4 = (1)* (line) + centerY 
    square = [x1',x2',y1',y2',x3',x4',y3',y4'];
    
end
%%%%%%%%%%%%%%%%
% Frame content specification for zooming in on the first quadrant of the
% plot, at time tau on a range of 0 to 1, demonstrating the simplified
% pan-zoom scheme. The extra argument "start_limits" allows this function
% to zoom from the end of the previous movie to the final position, rather
% than having its starting axis limits hard-coded
% function frame_info = animation_demo_panzoom(start_limits,frame_info, tau)
% 	
% 	% Declare ending axis limits
% 	end_limits = [0 .5*pi 0 1.1];
% 	
% 	% Set the axis limits to their corresponding values at fraction tau of
% 	% the way between their starting and ending values
% 	new_limits = panzoom(start_limits,end_limits,tau);
% 	axis(frame_info.ax,new_limits)
% 		
% 	% Declare a print method (in this case, print 150dpi png files of
% 	% frame_info.f, using the Painters renderer)
% 	frame_info.printmethod = @(dest) print(frame_info.f,'-dpng','-r 150','-painters',dest);
% 
end