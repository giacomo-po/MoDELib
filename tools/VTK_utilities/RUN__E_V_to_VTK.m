%% Takes data from the E and V directories and turns this into VTK data
%  Usage:
%  ======
%  1. set the variable path_to_E_and_V to a valid location where the E and
%     V directories are
%  2. define the variable "step" as a list of step numbers
%  3. define the output irectory for vtk data in the variable "vtk_outdir"
%
%  Info:
%  =====
%  - polygons are numbered from 1...
%  - vertex number are taken from the V file and +1 (Matlab one based)
%  - in ParaView use the "tube" filter for getting thicker, 3D-ish lines
%    (see the screenshot in "paraview_screenshot.png" for some information
%    about how to set up the visualization)
%  
%  ToDos:
%  ======
%  - replace 2-node polygons by only one polygon in the vtk file for
%    better performance
%  - add option for binary input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all


%% these are the variables to be user defined
path_to_E_and_V = 'ascii_data';
steps = [0 500 1000];
vtk_outdir = '.';     
vtk_basename = 'polygons';


% test if output dir exsits, create if it doesn't
if ~exist(vtk_outdir, 'dir')
    mkdir(vtk_outdir);
end

% format of Burgers vector string
bfmt = '%+06.3f %+06.3f %+06.3f';

%% get unique Burgers vector list of all simulation steps
idx = 1;
for step = steps
    fprintf('sorting Burgers vectors, step=%i\n',step)
    E=load(fullfile(path_to_E_and_V, ['/E/E_' num2str(step) '.txt']));
    for e=1:size(E,1)
        b=E(e,[3:5]);
        b_list{idx} = sprintf(bfmt,b);
        idx = idx +1;
    end
end
b_list = unique(b_list,'rows');

%% store data from E and V file into vertex and polygon cell array
for step = steps
    fprintf('writing VTK files, step=%i\n',step)
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
        vertex0 = round(E(e,1)) + 1;
        vertex1 = round(E(e,2)) + 1;
        polygon{e}.vertices = [vertex0, vertex1];
        polygon{e}.nvertices = 2;
        polygon{e}.b_index = find(ismember(b_list,sprintf(bfmt,b)));
    end
    
    vtk_file = sprintf('%03i_%s.vtk',step,vtk_basename);
    write_VTK_polygone_file( vtk_outdir, vtk_file , polygon, vertex);
end
