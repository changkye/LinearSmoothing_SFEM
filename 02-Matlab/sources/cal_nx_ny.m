function  [nx,ny] = cal_nx_ny(x,y,side)
	%--------------------------------------------------------
	%	compute the outside normal vector of sides of the element
	%	
	%	Synopsis:
	%	[nx,ny]=cal_nx_ny(x,y,side)
	%	
	%	variable description:
	%	nx - x direction normal of side
	%	ny - y direction normal of side
	%	side - the length of each side of the element
	%	x - x coord value of element
	%	y - y coord value of element
	%--------------------------------------------------------

	nsel = length(x);     %nsel n: number   s:side   el:element

	for ie = 1:nsel    
    	if ie == nsel
        	nx(nsel) =  (y(1)-y(nsel))/side(nsel); % x direction normal of nsel th side
        	ny(nsel) = -(x(1)-x(nsel))/side(nsel); % y direction normal of nsel th side
    	else
	        nx(ie) =  (y(ie+1)-y(ie))/side(ie);    % x direction normal of ie th side
    	    ny(ie) = -(x(ie+1)-x(ie))/side(ie);    % y direction normal of ie th side
    	end
	end

end