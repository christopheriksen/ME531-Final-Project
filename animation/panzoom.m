function intermediate_limits = panzoom(start_limits,end_limits,fraction)
% Find the axis limits that are FRACTION of the way between START_LIMITS
% and END_LIMITS, using a linear scale for centering and a log scale for
% zooming
%
%
% ========================
% Copyright Ross L. Hatton, 2011

	%%%%%%%%%
	% Limit centering

	% Get the center of each set of limits
	start_center = (start_limits(1:2:end)+start_limits(2:2:end))/2;
	end_center = (end_limits(1:2:end)+end_limits(2:2:end))/2;
	
	% take a weighted average
	intermediate_center = start_center*(1-fraction)+end_center*fraction;
	
	
	%%%%%%%%
	% Limit range
	
	% Get the range of each set of limits
	start_range = start_limits(2:2:end)-start_limits(1:2:end);
	end_range = end_limits(2:2:end)-end_limits(1:2:end);
	
	% take a logarithmic interpolation
	intermediate_range = exp(log(start_range)*(1-fraction) + log(end_range)*fraction);
	
    % If any of the limits are the same at beginning and end, anchor at
    % those points; otherwise use the intermediate center.
	% Build the new set of limits as the center +/- half the range
	
    % Placeholder
    intermediate_limits = zeros(size(start_limits));        
   	
	for i = 1:length(intermediate_center)
        
        % If the lower boundary of the ith dimension is the same, anchor
        % that boundary
        if start_limits(2*i-1) == end_limits(2*i-1)
            
            intermediate_limits(2*i-1) = start_limits(2*i-1);
            intermediate_limits(2*i) = start_limits(2*i-1)+intermediate_range(i);
            
        % If the upper boundary of the ith dimension is the same, anchor
        % that boundary
        elseif start_limits(2*i) == end_limits(2*i)
             
            intermediate_limits(2*i-1) = start_limits(2*i)-intermediate_range(i);
            intermediate_limits(2*i) = start_limits(2*i);
           
        % Otherwise, use the intermediate center
        else
            
            intermediate_limits(2*i-1) = intermediate_center(i)-.5*intermediate_range(i);
            intermediate_limits(2*i) = intermediate_center(i)+.5*intermediate_range(i);
        end
	end

end