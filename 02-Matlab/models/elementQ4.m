function [Nodes,Elements] = elementQ4(L,numEls)
	% 
	xn = numEls(1) + 1;
	yn = numEls(2) + 1;
	
	% nodal coordinates
	[x,y] = meshgrid(linspace(L(1,1),L(1,2),xn),linspace(L(2,1),L(2,2),yn));
	x = x'; y = y';
	Nodes = [x(:) y(:)];
	
	% element connectivities
	idx = reshape(1:prod([xn,yn]),[xn,yn]);
	v1 = idx(1:end-1,1:end-1);
	v2 = idx(1:end-1,2:end);
	v3 = idx(2:end,1:end-1);
	v4 = idx(2:end,2:end);
	v1 = v1(:); v2 = v2(:); v3 = v3(:); v4 = v4(:); 
	Elements = [v1 v3 v4 v2];
end