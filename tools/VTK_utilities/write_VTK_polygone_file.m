function write_VTK_polygone_file( vtk_dir, vtkfile, polygon, vertex )
%write_VTK_polygone_file Write a VTK file with dislocations as polygones

    
    Nv = size(vertex,2);  % total number of vertices
    Nd = size(polygon,2); % total number of dislocations
    
    % count total number of segments
    tmp = [polygon{:}];
    Nseg = sum([tmp.nvertices]);
    
    
    % tmp = cell2mat(polygon);
%     % the formatting string "burgers_format" is defined above
%     for nd=1:Nd
%         b = sprintf(burgers_format, polygon{nd}.b_vector);
%         idx = find(ismember(all_Burgers_vectors,b_vector));
%         sprintf('b=%s,  idx=%i,  all_Burgers_vectors=%s\n',b_vector,idx,all_Burgers_vectors{idx})
%     end

    % Header
    fid=fopen(fullfile(vtk_dir,vtkfile),'w');
    fprintf(fid,'%s\n','# vtk DataFile Version 1.0');
    fprintf(fid,'%s\n','Dislocation data');
    fprintf(fid,'%s\n','ASCII');
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','DATASET POLYDATA');
    fprintf(fid,'POINTS %i float\n',Nv);
    %
    for nv=1:Nv
        fprintf(fid, '%1.15e %1.15e %1.15e \n', vertex{nv}.coords);
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'LINES %i %i\n',Nd,Nd + Nseg);
    b_id = zeros(Nd,1); % ids of all possible burgers vectors for each polygon
    for p=1:Nd
        %fprintf('b_id=%i,  b=[%f, %f, %f]\n',b_id(p), b(1), b(2), b(3));
        
        fprintf(fid,'%i ',polygon{p}.nvertices);
        for nv=1:polygon{p}.nvertices
            fprintf(fid,'%i ', polygon{p}.vertex(nv)-1 );
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

