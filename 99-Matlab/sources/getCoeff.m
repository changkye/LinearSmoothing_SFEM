function [caa,cbb,caa1,ca1a,cbb1,cb1b] = getCoeff(va,vb,ell)
	% filename: getCoeff.m
	%
	% purpose: to compute the terms in the A matrix that relates the tensor
	% product of shape functions to the serendipity shape functions.
	%
	% Inputs:
	%       va, vb - (x,y) coordinates of the points a, b. Ordered
	%       non-repeating pair (a,b) forms the diagonal ends
	%
	%       a-1/b+1 ------------- b
	%           |                 |
	%           |                 |
	%           |                 |
	%           |                 |
	%           |                 |
	%           a ----------------a+1/b-1
	%
	% Note: in the above illustration, the points (a,b) - form the diagonal
	% pair
	%       ell - length of the diagonal
	%
	% Outputs:
	%       caa, cbb, c^{(a-1)a}, c^{a(a+1)}, c^{(b-1)b} and c^{b(b+1)}.
	%
	% Note: the output values are for the chosen diagonals and form the column
	% vector of the A matrix.
	%
	% ref: Quadratic serendipity finite elements on polygons using generalized
	% barycentric coordinates, Mathematics of Computation, v83, number290,
	% 2014, 2691--2716, Alexander Rand, Andrew Gillette and Chandrajit Bajaj
	%
	% Sundar, 2014
	%--------------------------------------------------------------------------
	
	% to compute the terms in the A matrix, either the x or the y coordinate
	% can be chosen. chosing either of them will yield the same result.
	% x coordinates of the point a
	va1a = va(1,1); vaa = va(2,1); vaa1 = va(3,1);

	% x coordinates of the point b
	vb1b = vb(1,1); vbb = vb(2,1); vbb1 = vb(3,1);

	% the intercepts, da and db. Ref: equation 5.2 and 5.3
	da = 1/ell*(va(1,1)*va(3,2) - va(3,1)*va(1,2))/(va(1,2)-va(3,2));

	db = 1/ell*(vb(3,1)*vb(1,2) - vb(1,1)*vb(3,2))/(vb(1,2)-vb(3,2));

	% s. Ref: equation 5.4
	s = 2/(2-(da+db)) ;

	% c^{aa} and c^{bb}. Ref: equation 5.9, 5.10 of the paper 
	caa = (-2-2*da)/(2-(da+db));
	cbb = (-2-2*db)/(2-(da+db));

	% c^{(a-1)a} and c^{a(a+1)}. Ref: equations 5.5 and 5.6
	% va1 - coordinates of the point (a-1), va2 - coordinates of the point a 
	% and va3 - coordinates of the point (a+1)
	Ate = [1 1;va1a vaa1];
	fte = [s; s*da*vaa];

	cte = Ate\fte;

	% c^{(a-1)a}
	ca1a = cte(1);

	%c^{a(a+1)}
	caa1 = cte(2);

	% c^{(b-1)b} and c^{b(b+1)}. Ref: equations 5.7 and 5.8
	% vb1 - coordinates of the point (b-1), vb2 - coordinates of the point b 
	% and vb3 - coordinates of the point (b+1)
	clear Ate fte cte
	Ate = [1 1;vb1b vbb1];
	fte = [s; s*db*vbb];

	cte = Ate\fte;

	% c^{(b-1)b}
	cb1b = cte(1);

	% c^{b(b+1)}
	cbb1 = cte(2);
end