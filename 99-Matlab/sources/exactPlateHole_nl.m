function [S,U] = exactPlateHole_nl(params,gpt,a)
	% Compute the exact stresses of the infinite plate with 
	% a circular centered hole problem
	%
	% Inputs:
	%   pt - point at which the solution needs to be computed
	%   a - radius of the hole...
	%--------------------------------------------------------------------------

	[theta,r] = cart2pol(gpt(1,1),gpt(1,2));

	c1t = cos(theta);
	c2t = cos(2*theta);
	c3t = cos(3*theta);
	c4t = cos(4*theta);
	s1t = sin(theta);
	s2t = sin(2*theta);
	s3t = sin(3*theta);
	s4t = sin(4*theta);
	fac1 = (a/r)^2;
	fac2 = fac1*fac1;

	S(1) = 1 - fac1*(1.5*c2t+c4t) + 1.5*fac2*c4t;
	S(2) =   - fac1*(0.5*c2t-c4t) - 1.5*fac2*c4t;
	S(3) =   - fac1*(0.5*s2t+s4t) + 1.5*fac2*s4t;

	mu   = params(1);	% shear modulus
	nu = (3*params(2)-2*params(1))/(2*(3*params(2)+params(1)));
	fac  = a/(8*mu);
	fac1 = 2*(a/r)^3;
	fac2 = 2*a/r;

	% plane strain
	k = 3 - 4*nu;
	% % plane stress
	% k = (3-nu)/(1+nu);
	
	U(1) = fac*(r/a*(k+1)*c1t + fac2*((k+1)*c1t + c3t) - fac1*c3t);
	U(2) = fac*(r/a*(k-3)*s1t + fac2*((1-k)*s1t + s3t) - fac1*s3t);
end