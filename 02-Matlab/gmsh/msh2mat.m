function msh2mat(bcType,probType)
	% gmsh msh to mat file
	format long; close all; clc;
    
	path(path,'./T3_to_T6');
    path(path,'../sources');

    elemType = 'T6';
    flag = 1;
    probPath = ['./' bcType '/' bcType '_' probType];
	in = fopen([probPath '.msh'],'r');
   
	tline = fgetl(in);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(in);
    end
    fclose(in);
    
    %----------------------------------------------------------------------------------------------
    % Get coordinates
    %----------------------------------------------------------------------------------------------
    % Find the tlines with texts on them
    models = tlines;
    txtLines = regexp(tlines,'.*Nodes','match','once');
    txtLineMask = ~cellfun(@isempty, txtLines);
    
    % Convert the non equation lines to the second numeric value
    for i = find(~txtLineMask)'
        Numbers = str2num(tlines{i});
        tlines{i} = Numbers;
    end
   
    % And the numbers in an array in the second cell
    txtLineNos = [find(txtLineMask); length(tlines)+1];
    for i = 1:length(txtLineNos)-1
        inds = txtLineNos(i)+1:txtLineNos(i+1)-1;
        nodes = tlines(inds);
        n=0;
        for j = 1:size(nodes,1)
            if size(nodes{j},2) == 3
                n = n + 1;
                blocks{n,1} = nodes{j};
            end
        end        
    end
    
    % get Nodes
    nnode = size(blocks,1);
    Nodes = zeros(nnode,2);
    for i = 1:nnode
        Nodes(i,:) = blocks{i}(1:end-1);
    end
    clear txtLines txtLineMask blocks inds nodes Numbers
    
    %----------------------------------------------------------------------------------------------
    % Get element connectivities
    %----------------------------------------------------------------------------------------------
    % Find the tlines with texts on them
    tlines = models;
    txtLines = regexp(tlines,'.*Elements','match','once');
    txtLineMask = ~cellfun(@isempty, txtLines);
    
    % Convert the non equation lines to the second numeric value
    for i = find(~txtLineMask)'
        Numbers = str2num(tlines{i});
        tlines{i} = Numbers;
    end
   
    % And the numbers in an array in the second cell
    txtLineNos = [find(txtLineMask); length(tlines)+1];
    for i = 1:length(txtLineNos)-1
        inds = txtLineNos(i)+1:txtLineNos(i+1)-1;
        nodes = tlines(inds);
        for j = 1:size(nodes,1)-2
            blocks{j,1} = nodes{j+2};
        end        
    end
    
    % get Elements
    nelem = size(blocks,1);
    Elements = zeros(nelem,3);
    for i = 1:nelem
        Elements(i,:) = blocks{i}(2:end);
    end
    
    %----------------------------------------------------------------------------------------------
    % Check
    %----------------------------------------------------------------------------------------------
    if flag == 1
        plotMesh(Nodes,Elements,'T3','k');
	   for in = 1:size(Nodes,1)
    	   xc = Nodes(in,1) + 0.001;
    	   yc = Nodes(in,2) - 0.003;
    	   text(xc,yc,num2str(in),'color','blue');
	   end
	   for iel = 1:size(Elements,1)
		  econ = Elements(iel,:);
		  nod = Nodes(econ,:);
		  text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
	   end
    end
    %----------------------------------------------------------------------------------------------
    % for T6 element
    %----------------------------------------------------------------------------------------------
    if strcmp(elemType,'T6')
        % create txt file for "triangulation_l2q.m"
        o_node = fopen([probPath '_nodes.txt'],'w');
        o_elem = fopen([probPath '_elements.txt'],'w');
        numNodes = size(Nodes,1); numEls = size(Elements,1);
        for i = 1:numNodes
            fprintf(o_node,'%f %f',Nodes(i,1),Nodes(i,2));
            if i~=numNodes
                fprintf(o_node,'\n');
            end
        end
        for i = 1:numEls
            fprintf(o_elem,'%d %d %d',Elements(i,1),Elements(i,2),Elements(i,3));
            if i~=numEls
                fprintf(o_elem,'\n');
            end
        end
        triangulation_l2q([probPath]);
        fclose(o_node); fclose(o_elem);
        
        clear o_node o_elem
        [Nodes,Elements] = tri_lin2quad([probPath]);
    end
    
    % save mat file
    save([probPath '.mat'],'Nodes','Elements');
    
    
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