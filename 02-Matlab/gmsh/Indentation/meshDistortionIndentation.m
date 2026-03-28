function meshDistortionIndentation(elemType,numEls)
	clc; close all;
	path(path,'./sources');

	L = [0 2; 0 1];
	mesh = str2func(['element' elemType]);
	[Nodes,Elements] = mesh(L,numEls);
	btmNodes = find(Nodes(:,2) == min(Nodes(:,2)));
	topNodes = find(Nodes(:,2) == max(Nodes(:,2)));
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	% mesh distortion
	air = 0.4;
	dx = L(1,2)/numEls(1); dy = L(2,2)/numEls(2);
	xn = numEls(1) + 1; yn = numEls(2) + 1;
	for i = 1:size(Nodes,1)
		r = random('beta',1,1);
		r = air*(2*r - 1);
		if sum(i == unique([btmNodes',topNodes',leftNodes',rightNodes'])) == 0
			Nodes(i,1) = Nodes(i,1) + dx*r;
			Nodes(i,2) = Nodes(i,2) + dy*r;
		end
	end
	flag = 1;
	if flag == 1
		clf; axis equal; axis on; hold on;
    	plotMesh(Nodes,Elements,elemType,'k');
% 		for in = 1:size(Nodes,1)
%     		xc = Nodes(in,1) + 0.001;
%     		yc = Nodes(in,2) - 0.003;
%     		text(xc,yc,num2str(in),'color','blue');
% 		end
% 		for iel = 1:size(Elements,1)
% 			econ = Elements(iel,:);
% 			nod = Nodes(econ,:);
% 			text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
% 		end
	end
	
	save(['./gmsh/Indentation/Indentation_' num2str(numEls(1)) 'x' num2str(numEls(2)) '.mat'],...
		'Nodes','Elements');

end