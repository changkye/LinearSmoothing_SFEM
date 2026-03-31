function [phi,dphi] = wachspress2D(v,x)
	% Inputs
	% v = [x1 y1; x2 y2...] - the n vertices of the polygon in ccw direction
	%  x = [x(1) x(2)] - point at which shape function is required...
	% 
	% Outputs
	% phi - output basis functions
	% dphi - output gradient of shape functions
	% 
	% Floater, Gillette and Sukumar 2013
	% -----------------------------------------------------------------------

	n = size(v,1); w = zeros(n,1);
	R = zeros(n,2);
	phi = zeros(n,1);
	dphi = zeros(n,2);

	% get the normals...
	un = zeros(n,2);
	for i = 1:n
		d = v(mod(i,n)+1,:) - v(i,:);
		un(i,:) = [d(2) -d(1)]/norm(d);
	end

	p = zeros(n,2);
	for i = 1:n
		h = dot(v(i,:)-x,un(i,:));
		p(i,:) = un(i,:)/h;
	end

	for i = 1:n
		im1 = mod(i-2,n) + 1;
		w(i) = det([p(im1,:); p(i,:)]);
		R(i,:) = p(im1,:) + p(i,:);
	end

	wsum = sum(w);
	phi = w/wsum;

	phiR = phi'*R;
	for k = 1:2
		dphi(:,k) = phi.*(R(:,k)-phiR(:,k));
	end
end
