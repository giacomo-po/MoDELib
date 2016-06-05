path_to_E_and_V = 'refdata';
step = 1;
vtk_dir = '.';     
vtk_file = sprintf('%03i_polygones.vtk',step);

% %% get Burgers vectors of all simulation steps
% Burgers=[];
% for step = [1]
%     % read vertex data
%     V=load(fullfile(path_to_E_and_V, ['/V/V_' num2str(step) '.txt']);
%     
%     % read edge data
%     E=load(fullfile(path_to_E_and_V, ['/E/E_' num2str(step) '.txt']);
%    
%     for e=1:size(E,1)
%         b=E(e,[3:5]);
%         Burgers=[Burgers; b];
%     end
% end
% Burgers = unique(Burgers,'rows');


vertex{1}.coords = [ 10., 10., 0. ];
vertex{2}.coords = [ 25., 10., 0. ];
polygon{1}.vertex = [1, 2];
polygon{1}.nvertices = 2;
polygon{1}.b_index = 1;

vertex{3}.coords = [ 10., 15., 0. ];
vertex{4}.coords = [ 20., 18., 8. ];
vertex{5}.coords = [ 25., 20., 10. ];
polygon{2}.vertex = [3, 4, 5];
polygon{2}.nvertices = 3;
polygon{2}.b_index = 2;

write_VTK_polygone_file( vtk_dir, vtk_file , polygon, vertex);
