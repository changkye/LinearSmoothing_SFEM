function [kk,bb] = feaplyc2_nl(model,uu,scale)
	%----------------------------------------------------------
	%  Purpose:
	%     Apply constraints to matrix equation [kk]{x}={ff}
	%
	%  Synopsis:
	%     [kk,ff]=feaplybc(kk,ff,bcdof,bcval)
	%
	%  Variable Description:
	%     kk - system matrix before applying constraints 
	%     ff - system vector before applying constraints
	%     bcdof - a vector containging constrained d.o.f
	%     bcval - a vector containing contained value 
	%     For example, there are constraints at clamp edge and 
	%     their constrained values are 0.0  
	%-----------------------------------------------------------
	kk = full(model.K); rr = full(model.R); ff = model.F; 
	bcdof = model.bcdof; bcval = model.bcval;
 	    
 	n = length(bcdof);
 	sdof = size(kk,1);
    if length(bcval) ~= n
        error('feaplyc2_nl:BCSizeMismatch', ...
            'Length of bcval (%d) must match length of bcdof (%d).', length(bcval), n);
    end

 	bb = scale*ff - rr;
	for i = 1:n
    	c = bcdof(i);
    	kk(c,1:sdof) = 0.;
    	kk(c,c) = 1.;
%         kk(c,c) = 1e7;
    	bb(c) = scale*bcval(i) - uu(c);
 	end

 end
