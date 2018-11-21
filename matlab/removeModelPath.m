function removeModelPath(dirName)

fList = dir(dirName);

for file = {fList.name}
    
    fileName=file{1};
    
    
    if(length(fileName)>2 || fileName=="IO")
        
        if isfolder([dirName '/' fileName])
            removeModelPath([dirName '/' fileName])
        end
        
        if(fileName(end-1:end)=='.h')
            fid = fopen([dirName '/' fileName],'r');
            while true
                line=fgetl(fid);
                if ~ischar(line)
                    break
                end
                
                foundIncludeMode=contains(line,"#include <model");
                foundDotH=contains(line,".h>");
                if  foundIncludeMode && foundDotH
                
                    foundLess=strfind(line,'<');
                    foundGreater=strfind(line,'>');
                    oldPath=line(foundLess+1:foundGreater-1);
                    if contains(oldPath,'/')
                    lastSlash=strfind(oldPath,'/');
                    newPath=oldPath(lastSlash(end)+1:end);
                    [oldPath ' ---> ' newPath]
                    command=['sed -i ' '''' '''' ' ' '''' 's,'  oldPath   ',' newPath ',' '''' ' ' dirName '/' fileName ]
                    system(command)
                    end
                end
                
            end
            fclose(fid);
        end
        
    end
    
    
    
end