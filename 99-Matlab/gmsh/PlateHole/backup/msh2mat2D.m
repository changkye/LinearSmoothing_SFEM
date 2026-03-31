function model = msh2mat2D(model,probType)
	% gmsh msh to mat file
	% 
    % clc;
    % restoredefaultpath;
	path(path,'./T3_to_T6');
    
	in = fopen(['./' model '/' model '_' probType '.msh'],'r');
   
	cellarray = textscan(in,'%s');
	cellno = 7;

	etc1 = str2num(cellarray{1}{cellno});
	cellno = cellno + 1;
	etc2 = str2num(cellarray{1}{cellno});
	cellno = cellno + 1;
	etc3 = str2num(cellarray{1}{cellno});
	cellno = cellno + 2;

    dummy1 = zeros(etc1,4); ndummy1 = size(dummy1,2);
    dummy2 = zeros(etc2,10); ndummy2 = size(dummy2,2);
    dummy3 = zeros(etc3,14); ndummy3 = size(dummy3,2);
	for i = 1:etc1
		for j = 1:ndummy1
        	dummy1(i,j) = str2num(cellarray{1}{cellno});
            cellno = cellno + 1;
        end
        cellno = cellno + 1;
	end
	cellno = cellno + 1;
	for i = 1:etc2
		for j = 1:ndummy2
			dummy2(i,j) = str2num(cellarray{1}{cellno});
			cellno = cellno + 1;
		end
		cellno = cellno + 1;
	end
% 	cellno = cellno + 1;
	for i = 1:etc3
		for j = 1:ndummy3
			dummy3(i,j) = str2num(cellarray{1}{cellno});
			cellno = cellno + 1;
		end
		cellno = cellno + 1;
	end
	
	% coordinates
	cellno = cellno + 1;
	dummy = str2num(cellarray{1}{cellno});
	cellno = cellno + 1;
	numNode = str2num(cellarray{1}{cellno});
	xy = zeros(numNode,4);
    cellno = cellno + 2;
    
	for i = 1:numNode+dummy
		for j = 1:4
            cellno = cellno + 1;
			xy(i,j) = str2num(cellarray{1}{cellno});
        end
    end
    Nodes = zeros(numNode,2);
    n = 0;
   	for i = 1:numNode+dummy
   		if find(xy(i,end)==0)
   			n = n + 1;
   			Nodes(n,:) = xy(i,2:end-1);
   		end
    end
    
    % connectivity
    cellno = cellno + 4;
    numEls = str2num(cellarray{1}{cellno});
    Elements = zeros(numEls,3);
    cellno = cellno + 5;
    for i = 1:numEls
        for j = 1:3
        	cellno = cellno + 1;
        	Elements(i,j) = str2num(cellarray{1}{cellno});
        end
        cellno = cellno + 1;
    end
    Elements = Elements - 1;
    fclose(in);
        
    if strcmp(model.elemType,'T6')
        % create txt file for "triangulation_l2q.m"
    	o_node = fopen(['./gmsh/' model.bcType '/' ...
            model.bcType '_' probType '_nodes.txt'],'w');
    	o_elem = fopen(['./gmsh/' model.bcType '/' ...
            model.bcType '_' probType '_elements.txt'],'w');
    	numNodes = size(Nodes,1); numEls = size(Elements,1);
    	for i = 1:numNodes
    		fprintf(o_node,'%f  %f',Nodes(i,1),Nodes(i,2));
    		if i~=numNodes
    			fprintf(o_node,'\n');
    		end
    	end
    	for i = 1:numEls
    		fprintf(o_elem,'%d  %d  %d',Elements(i,1),Elements(i,2),Elements(i,3));
    		if i~=numEls
    			fprintf(o_elem,'\n');
    		end
        end
		triangulation_l2q(['./gmsh/' model.bcType '/' model.bcType '_' probType]);
		fclose(o_node); fclose(o_elem);
		% save mat file
		clear o_node o_elem
		[Nodes,Elements] = tri_lin2quad(['./gmsh/' model.bcType '/' model.bcType '_' probType]);
    end
    
    % save mat file
    % save([infile '.mat'],'Nodes','Elements');
    model.Nodes = Nodes;
    model.Elements = Elements;
    
end 
% 
function [Nodes,Elements] = tri_lin2quad(prefix)
	node_l2q_filename = strcat(prefix,'_l2q_nodes.txt');
    element_l2q_filename = strcat(prefix,'_l2q_elements.txt');
    
    % read coordinates
    [dim_num,node_num1] = r8mat_header_read(node_l2q_filename);
    node_xy1(1:dim_num,1:node_num1) = r8mat_data_read(node_l2q_filename,dim_num,node_num1);
    
    % read connectivities    
    [triangle_order1,triangle_num] = i4mat_header_read(element_l2q_filename);
    triangle_node1(1:triangle_order1,1:triangle_num) = i4mat_data_read(...
        element_l2q_filename,triangle_order1,triangle_num);
    
    Nodes = node_xy1';
    Elements = triangle_node1';
end