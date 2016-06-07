function write_VTK_polygone_file( vtk_dir, vtkfile, polygon, vertex )
%write_VTK_polygone_file Write a VTK file with dislocations as polygones
   
    Nv = size(vertex,2);  % total number of vertices
    Nd = size(polygon,2); % total number of dislocations
    
    % count total number of segments
    tmp = [polygon{:}];
    Nseg = sum([tmp.nvertices]);
    
     fid=fopen(fullfile(vtk_dir,vtkfile),'w');
    fprintf(fid,'%s\n','# vtk DataFile Version 1.0');
    fprintf(fid,'%s\n','Dislocation data');
    fprintf(fid,'%s\n','ASCII');
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','DATASET POLYDATA');
    fprintf(fid,'POINTS %i float\n',Nv);
    
    % find valid vertex coordinates for filling up empty nodeIDs
    for nv=1:Nv
        if isempty(vertex{nv}); 
            continue;
        end
        dummy_coords = vertex{nv}.coords;
        break;
    end
    
    for nv=1:Nv
        if isempty(vertex{nv}); 
            coords = dummy_coords; 
        else
            coords = vertex{nv}.coords;
        end
        fprintf(fid, '%1.8e %1.8e %1.8e \n', coords);
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'LINES %i %i\n',Nd,Nd + Nseg);
    for p=1:Nd
        fprintf(fid,'%i ',polygon{p}.nvertices);
        for nv=1:polygon{p}.nvertices
            fprintf(fid,'%i ', polygon{p}.vertices(nv)-1 );
        end
        fprintf(fid,'\n');
    end

    fprintf(fid,'\nCELL_DATA %i\n',Nd);
    fprintf(fid,'SCALARS Burgers_vector_ID int 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    for p=1:Nd
        fprintf(fid, '%i \n', polygon{p}.b_index);
    end
    
    fprintf(fid,'SCALARS numbers int 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    for p=1:Nd
        fprintf(fid, '%i \n', p);
    end
    
    fclose(fid);
end

