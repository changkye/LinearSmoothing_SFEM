function [Nodes,Elements] = elementQ9(L,numEls)
	% 
	xn = 2*numEls(1) + 1;
	yn = 2*numEls(2) + 1;
	
	% nodal coordinates
	[x,y] = meshgrid(linspace(L(1,1),L(1,2),xn),linspace(L(2,1),L(2,2),yn));
	x = x'; y = y';
	Nodes = [x(:) y(:)];
	
	% element connectivity
	ic = 2; jc = 2*xn;
	idx = [1 3 2*xn+3 2*xn+1 2 xn+3 2*xn+2 xn+1 xn+2];
	Elements = makeConnectivity(model,ic,jc,idx);

end
%%
function element = makeConnectivity(numEls,ic,jc,idx)
	% 
    num_u = numEls(1); num_v = numEls(2);
    
	if nargin < 4
        disp(['Not enough parameters specified for make_elem function'])
    end

    inc = zeros(1,size(idx,2));
    e = 1;
    element = zeros(num_u*num_v,size(idx,2));

    for row = 1:num_v
        for col = 1:num_u
            element(e,:) = idx + inc;
            inc = inc + ic;
            e = e + 1;
        end
        inc = row*jc;
    end
end