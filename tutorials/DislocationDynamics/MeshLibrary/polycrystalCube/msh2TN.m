clear all
close all
clc

addpath('../../../../matlab');

filename='cube50.msh';
meshID=3;
centerAndScale=1;
scaleFactor=4000;

fid = fopen(filename);
tline = fgetl(fid);
while ischar(tline)
    
    if strfind(tline,'$Nodes')
        tline = fgetl(fid);
        nNodes=eval(tline)
        Nodes=zeros(nNodes,4);
        for n=1:nNodes
            tline = fgetl(fid);
            row=eval(['[' tline ']']);
            Nodes(n,:)=row;
        end
    end
    
    
    if strfind(tline,'$Elements')
        tline = fgetl(fid);
        nElements=eval(tline)
        Tets=[];
        for n=1:nElements
            tline = fgetl(fid);
            row=eval(['[' tline ']']);
            if(row(2)==4)
                region1=row(4);
                region2=row(5);
                if(region1~=region2)
                    error('region mismatch')
                end
                Tets=[Tets;[row(1) row(7:10) region1]];
            end
        end
    end
        
    tline = fgetl(fid);
end
fclose(fid);

%%
if centerAndScale
    xyz=Nodes(:,2:4);
    xyz=xyz-repmat(mean([min(xyz);max(xyz)]),size(xyz,1),1);
    xyz=xyz*scaleFactor;
    Nodes(:,2:4)=xyz;
end

%return
%% Write N_meshID.txt
nodeFile = fopen(['N_' num2str(meshID) '.txt'],'w');
nodeFormat='%i %.15e %.15e %.15e \n';
for n=1:size(Nodes,1)
fprintf(nodeFile,nodeFormat,Nodes(n,:));
end
fclose(nodeFile);

tetFile = fopen(['T_' num2str(meshID) '.txt'],'w');
tetFormat='%i %i %i %i %i %i \n\n';
for n=1:size(Tets,1)
fprintf(tetFile,tetFormat,Tets(n,:));
end
fclose(tetFile);

%%
grainIDs=unique(Tets(:,end));
polyFile = fopen('polyCrystal.txt','w');
fprintf(polyFile,'materialFile=../../MaterialsLibrary/Cu.txt;\n');
polyFormat='%.15e %.15e %.15e';
for g=1:length(grainIDs)
    phi=rand(1)*2*pi;
    theta=rand(1)*pi;
    a=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]';
    R=angleAxis(a,rand(1)*2*pi);
    %R=eye(3);
fprintf(polyFile,['C2G' num2str(grainIDs(g)) '=']);
fprintf(polyFile,[polyFormat '\n'],R(1,:));
fprintf(polyFile,[polyFormat '\n'],R(2,:));
fprintf(polyFile,[polyFormat ';\n'],R(3,:));
%fprintf(polyFile,';\n');
end
fclose(polyFile);

