function triangulation_l2q(prefix)
    %*****************************************************************************80
    %
    % MAIN is the main program for TRIANGULATION_L2Q.
    %
    %  Discussion:
    %
    %    TRIANGULATION_L2Q makes a quadratic triangulation from a linear one.
    %
    %    Thanks to Zhu Wang for pointing out a problem caused by a change
    %    in the ordering of elements in the triangle neighbor array, 25 August 2010.
    %
    %  Usage:
    %
    %    triangulation_l2q ( 'prefix' )
    %
    %    where 'prefix' is the common filename prefix:
    %
    %    * 'prefix'_nodes.txt contains the node coordinates,
    %    * 'prefix'_elements.txt contains the element definitions.
    %    * 'prefix'_l2q_nodes.txt will contain the quadratic node coordinates,
    %    * 'prefix'_l2q_elements.txt will contain the quadratic element definitions.
    %
    %  Licensing:
    %
    %    This code is distributed under the GNU LGPL license.
    %
    %  Modified:
    %
    %    25 August 2010
    %
    %  Author:
    %
    %    John Burkardt
    %
    timestamp();

    fprintf(1,'\n');
    fprintf(1,'TRIANGULATION_L2Q\n');
    fprintf(1,'  MATLAB version\n');
    fprintf(1,'  Read a "linear" triangulation and\n');
    fprintf(1,'  write out a "quadratic" triangulation.\n');
    fprintf(1,'\n');
    fprintf(1,'  Read a dataset of NODE_NUM1 points in 2 dimensions.\n');
    fprintf(1,'  Read an associated triangulation dataset of TRIANGLE_NUM \n');
    fprintf(1,'  triangles which uses 3 nodes per triangle.\n');
    fprintf(1,'\n');
    fprintf(1,'  Create new nodes which are triangle midpoints,\n');
    fprintf(1,'  generate new node and triangulation data for\n');
    fprintf(1,'  quadratic 6-node triangles, and write them out.\n');
    %
    %  The command line argument is the common filename prefix.
    %
    if nargin < 1

        fprintf(1,'\n');
        fprintf(1,'TRIANGULATION_L2Q:\n');

        prefix = input('Please enter the filename prefix:');

    end
    %
    %  Create the filenames.
    %
    node_filename = strcat(prefix,'_nodes.txt');
    element_filename = strcat(prefix,'_elements.txt');
    node_l2q_filename = strcat(prefix,'_l2q_nodes.txt');
    element_l2q_filename = strcat(prefix,'_l2q_elements.txt');
    %
    %  Read the data.
    %
    [dim_num,node_num1] = r8mat_header_read(node_filename);

    fprintf(1,'\n');
    fprintf(1,'  Read the header of "%s".\n',node_filename);
    fprintf(1,'\n');
    fprintf(1,'  Spatial dimension DIM_NUM = %d\n',dim_num);
    fprintf(1,'  Number of points NODE_NUM1  = %d\n',node_num1);

    node_xy1(1:dim_num,1:node_num1) = r8mat_data_read(node_filename,dim_num,node_num1);

    fprintf(1,'\n');
    fprintf(1,'  Read the data in "%s".\n', node_filename);

    r8mat_transpose_print_some(dim_num,node_num1,node_xy1,1,1,5,5,...
        '  5 by 5 portion of data read from file:');
    %
    %  Read the triangulation data.
    %
    [triangle_order1,triangle_num] = i4mat_header_read(element_filename);

    if triangle_order1 ~= 3
        fprintf(1,'\n');
        fprintf(1,'TRIANGULATION_L2Q - Fatal error!\n');
        fprintf(1,'  Data is not for a 3-node triangulation.\n');
        error('TRIANGULATION_L2Q - Fatal error!');
    end

    fprintf(1,'\n');
    fprintf(1,'  Read the header of ""%s".\n', element_filename);
    fprintf(1,'\n');
    fprintf(1,'  Triangle order = %d\n', triangle_order1);
    fprintf(1,'  Number of triangles TRIANGLE_NUM  = %d\n', triangle_num);

    triangle_node1(1:triangle_order1,1:triangle_num) = i4mat_data_read(...
        element_filename,triangle_order1,triangle_num);

    fprintf(1,'\n' );
    fprintf(1,'  Read the data in ""%s".\n', element_filename);

    i4mat_transpose_print_some(triangle_order1,triangle_num,triangle_node1,...
        1,1,triangle_order1,10,'  3 by 10 portion TRIANGLE_NODE1:');
    %
    %  Detect and correct 0-based indexing.
    %
    triangle_node1 = mesh_base_one(node_num1,triangle_order1,triangle_num,triangle_node1);
    %
    %  Determine the number of midside nodes that will be added.
    %
    boundary_num = triangulation_order3_boundary_edge_count(triangle_num,triangle_node1);

    interior_num = (3*triangle_num-boundary_num)/2;
    edge_num = interior_num + boundary_num;
    fprintf(1,'\n');
    fprintf(1,'  Number of midside nodes to add = %d\n',edge_num);
    %
    %  Build the triangle neighbor array.
    %
    triangle_neighbor = triangulation_neighbor_triangles(triangle_order1,...
        triangle_num,triangle_node1);

    i4mat_transpose_print(3,triangle_num,triangle_neighbor,'  Triangle_neighbor');
    %
    %  Create the midside nodes.
    %
    triangle_order2 = 6;
    triangle_node2(1:3,1:triangle_num) = triangle_node1(1:3,1:triangle_num);
    triangle_node2(4:6,1:triangle_num) = -1;
    node_xy2(1:2,1:node_num1) = node_xy1(1:2,1:node_num1);

    node_num2 = node_num1;

    fprintf(1,'\n');
    fprintf(1,'  Generate midside nodes\n');
    fprintf(1,'\n');

    for triangle = 1:triangle_num

        for i = 1:3
            %
            %  CORRECTION #1 because element neighbor definition changed.
            %
            iii = i4_wrap(i+2,1,3);
            triangle2 = triangle_neighbor(iii,triangle);
    
            if (0 < triangle2 && triangle2 < triangle)
                continue
            end
    
            ip1 = i4_wrap(i+1,1,3);
    
            k1 = triangle_node2(i,triangle);
            k2 = triangle_node2(ip1,triangle);
    
            node_num2 = node_num2 + 1;
            node_xy2(1:dim_num,node_num2) = 0.5*(node_xy1(1:dim_num,k1)+node_xy1(1:dim_num,k2));
    
            fprintf(1,'%8d  %14f  %14f\n',node_num2,node_xy2(1:dim_num,node_num2));
    
            triangle_node2(3+i,triangle) = node_num2;
    
            if 0 < triangle2
                for ii = 1:3
                    %
                    %  CORRECTION #2 because element neighbor definition changed.
                    %
                    iii = i4_wrap(ii+2,1,3);
                    if triangle_neighbor(iii,triangle2) == triangle
                        triangle_node2(ii+3,triangle2) = node_num2;
                    end
                end
            end
        end
    end

    i4mat_transpose_print(triangle_order2,triangle_num,triangle_node2, ...
        '  TRIANGLE_NODE2' );
    %
    %  Write out the node and triangle data for the quadratic mesh.
    %
    r8mat_transpose_print(dim_num,node_num2,node_xy2,'  NODE_XY2:' );

    r8mat_write(node_l2q_filename,dim_num,node_num2,node_xy2);

    i4mat_write(element_l2q_filename,triangle_order2,triangle_num,...
        triangle_node2);
    %
    %  Terminate.
    %
    fprintf(1,'\n' );
    fprintf(1,'TRIANGULATION_L2Q\n');
    fprintf(1,'  Normal end of execution.\n');

    fprintf(1,'\n');
    timestamp();

    return
end