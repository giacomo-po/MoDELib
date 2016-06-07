% polygons are numbered from 1...
% vertex number are taken from the V file and +1 (Matlab one based)

%% TODO: replace 2-node polygons by only one polygon
clear all

path_to_E_and_V = '.';
step = 1;
vtk_dir = '.';     

%% get Burgers vectors of all simulation steps
Burgers=[];
for step = [0:10:100]
    fprintf('step=%i\n',step)
    % read vertex data
    V=load(fullfile(path_to_E_and_V, ['/V/V_' num2str(step) '.txt']));
    
    % read edge data
    E=load(fullfile(path_to_E_and_V, ['/E/E_' num2str(step) '.txt']));
   
    for v=1:size(V,1)
        nodeID = V(v,1);
        coords = V(v,2:4);
        vertex{nodeID+1}.coords = [coords(1), coords(2), coords(3)]; 
    end
    
    for e=1:size(E,1)
        b=E(e,[3:5]);
        bstring = sprintf('%1.1f %1.1f %1.1f',b);
        Burgers=[Burgers; b];
        vertex0 = round(E(e,1)) + 1;
        vertex1 = round(E(e,2)) + 1;
        polygon{e}.vertices = [vertex0, vertex1];
        polygon{e}.nvertices = 2;
        polygon{e}.b_index = 1;
    end
    
    vtk_file = sprintf('%03i_polygones.vtk',step);
    write_VTK_polygone_file( vtk_dir, vtk_file , polygon, vertex);

end
Burgers = unique(Burgers,'rows');
% 
% find(norm(Burgers-b)<1e-6)
