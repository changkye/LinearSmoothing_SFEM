function [Nodes,Elements] = elementQ8(L,numEls)
	% 
	xn = numEls(1) + 1;
	yn = numEls(2) + 1;
	
	% nodal coordinates
	for j = 1:yn
		for i = 1:xn
			elem = (3*numEls(1)+2)*(j-1)+(2*i-1);
			Nodes(elem,1) = (L(1,2)/numEls(1))*(i-1);
        	Nodes(elem,2) = (L(2,2)/numEls(2))*(j-1);
		end
	end
	
	% element connectivity
	index = 1;
	for i = 1:numEls(2)
		for k = 1:numEls(1)
        	Elements(index,1) = (3*numEls(1)+2)*(i-1) + (2*k-1);
        	Elements(index,2) = (3*numEls(1)+2)*(i-1) + (2*k+1);
        	Elements(index,3) = (3*numEls(1)+2)*i+(2*k+1);
        	Elements(index,4) = (3*numEls(1)+2)*i+(2*k-1);
        	Elements(index,5) = mean(Elements(index,[1,2])); 
        	Elements(index,6) = (3*numEls(1)+2)*i - numEls(1) + k;
        	Elements(index,7) = mean(Elements(index,[3,4])); 
        	Elements(index,8) = (3*numEls(1)+2)*i - (numEls(1)+1)+ k;
        	index = index + 1;
		end
	end

end