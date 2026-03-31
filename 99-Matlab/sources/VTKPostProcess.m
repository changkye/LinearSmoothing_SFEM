function VTKPostProcess(vtu,model,result)
    % Copyright (C) 2007 Garth N. Wells
    %
    % Write VTK post-processing files
    %       
    % Modified by VP Nguyen for the IGA FEM code
    %
    % Modified by Changkye Lee, School of Engineering,
    % iMAM, Cardiff University,
    % LeeC15@cardiff.ac.uk.
    
    DIM = size(model.Nodes,2);

    % Output files
    out = fopen([vtu],'w');

    if strcmp(model.elemType,'Q4')
        numVertexesPerCell = 4;
        VTKCellCode = 9;
    elseif strcmp(model.elemType,'T3') | strcmp(model.elemType,'T4') 
        numVertexesPerCell = 3;
        VTKCellCode = 5;
    elseif strcmp(model.elemType,'T6')
        numVertexesPerCell = 6;
        VTKCellCode = 22;
    elseif strcmp(model.elemType,'B8')
        numVertexesPerCell = 8;
        VTKCellCode = 12;
    elseif strcmp(model.elemType,'H4') | strcmp(model.elemType,'H4b')
        numVertexesPerCell = 4;
        VTKCellCode = 10;
    else
        error('Element type not known (VTKPostProcess)')
    end
    Nodes = model.Nodes;
    Elements = model.Elements;
    numNodes = size(Nodes,1);
    numEls = size(Elements,1);
    
    dof_per_vertex = 3;

    % Write headers
    fprintf(out,'<VTKFile type="UnstructuredGrid" version="0.1">\n');
    fprintf(out,'\t<UnstructuredGrid>\n');
    fprintf(out,'\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">\n',...
        numNodes,size(model.Elements,1));

    % Write nodal coordinates
    fprintf(out,'\t\t\t<Points>\n');
    fprintf(out,'\t\t\t\t<DataArray type="Float64" NumberOfComponents="%d" format="ascii">\n',...
        dof_per_vertex);
    if DIM == 2
        Nodes = cat(2,Nodes,repelem(0,size(Nodes,1),1));
    end
    
    for i = 1:numNodes
        fprintf(out,'\t\t\t\t\t');
        fprintf(out,'%f\t',Nodes(i,:));
        fprintf(out,'\n');
    end
    fprintf(out,'\t\t\t\t</DataArray>\n');
    fprintf(out,'\t\t\t</Points>\n');

    % Write Element connectivities
    fprintf(out,'\t\t\t<Cells>\n');
    % Print cell connectivity
    fprintf(out,'\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for i = 1:numEls
        fprintf(out,'\t\t\t\t\t');
        fprintf(out,'%d\t',Elements(i,1:numVertexesPerCell)-1);
        fprintf(out,'\n');
    end 
    fprintf(out,'\t\t\t\t</DataArray>\n');

    % Write offsets
    fprintf(out,'\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n');
    offset = 0;
    for i = 1:numEls
        offset = offset + numVertexesPerCell;
        fprintf(out,'\t\t\t\t\t%d\n',offset);
    end 
    fprintf(out,'\t\t\t\t</DataArray>\n');

    % Write cell code
    fprintf(out,'\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n');
    for i = 1:numEls
        fprintf(out,'\t\t\t\t\t%d\n',VTKCellCode);
    end 
    fprintf(out,'\t\t\t\t</DataArray>\n');
    fprintf(out,'\t\t\t</Cells>\n');

    % Write Displacement
    fprintf(out,'\t\t\t<PointData Vectors="U">\n');
    fprintf(out,'\t\t\t\t<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii"> \n');
    U = zeros(numNodes,3);
    if DIM == 2
        U = [result.uu(1:2:end), result.uu(2:2:end), zeros(numNodes,1)];
    else
        U = [result.uu(1:3:end), result.uu(2:3:end), result.uu(3:3:end)];
    end
    for i = 1:numNodes
        fprintf(out,'\t\t\t\t\t');
        fprintf(out,'%f\t',U(i,:));
        fprintf(out,'\n');
    end 
    fprintf(out,'\t\t\t\t</DataArray>\n');
    
    % Write strain
    if isfield(result,'E') == 1
        fprintf(out,'\t\t\t\t<DataArray type="Float64" Name="E" NumberOfComponents="%d" format="ascii">\n',size(result.E,2));
        for i = 1:numNodes
            fprintf(out,'\t\t\t\t\t');
            fprintf(out,'%f\t',result.E(wkInd(i),:));
            fprintf(out,'\n');
        end
        fprintf(out,'\t\t\t\t</DataArray>\n');
    end

    % Write stress
    if isfield(result,'S') == 1
        fprintf(out,'\t\t\t\t<DataArray type="Float64" Name="S" NumberOfComponents="%d" format="ascii">\n',size(result.S,2));
        for i = 1:numNodes
            fprintf(out,'\t\t\t\t\t');
            fprintf(out,'%f\t',result.S(wkInd(i),:));
            fprintf(out,'\n');
        end
        fprintf(out,'\t\t\t\t</DataArray>\n');
    end

    % Write von Mises stress
    if isfield(result,'vM') == 1
        fprintf(out,'\t\t\t\t<DataArray type="Float64" Name="von Mises" NumberOfComponents="1" format="ascii">\n');
        for i = 1:numNodes
            fprintf(out,'\t\t\t\t\t%f\n',result.vM(wkInd(i)));   
        end
        fprintf(out,'\t\t\t\t</DataArray>\n');
    end

    % Finishing
    fprintf(out,'\t\t\t</PointData>\n');
    fprintf(out,'\t\t</Piece>\n');
    fprintf(out,'\t</UnstructuredGrid>\n');
    fprintf(out,'</VTKFile>\n');
    fclose(out);
    
end
