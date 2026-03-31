function [side] = cal_side(x,y)
	%--------------------------------------------------------
	%	purpose:
	%	compute the length of each side of the element
	%	
	%	Synopsis:
	%	[side]=cal_side(x,y)
	%	
	%	variable description:
	%	side - the length of each side of the element
	%	x - x coord value of element
	%	y - y coord value of element
	%--------------------------------------------------------

	nsel = length(x);     %nsel n: number   s:side   el:element

	for ie = 1:nsel    
    	if ie == nsel
    		% length of side nsel        
	        side(nsel) = sqrt((x(nsel)-x(1))^2 + (y(nsel)-y(1))^2);        
    	else
			% length of side ie
        	side(ie) = sqrt((x(ie)-x(ie+1))^2 + (y(ie)-y(ie+1))^2);
    	end
	end

end