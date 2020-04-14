function printMatrixToFile(fid,M,label)
fprintf(fid,[label '=']);

format='';
for(c=1:size(M,2))
    format=[format '%1.15e '];
end

for(k=1:size(M,1))
    if k<size(M,1)
        fprintf(fid,[format '\n'],M(k,:));
    else
        fprintf(fid,[format ';\n'],M(k,:));
    end
end
end