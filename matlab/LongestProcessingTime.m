clear
clc

m=10; % number of processors
wM=zeros(m,1); % initial weights of the processors

% add a few heavy particles
np=7; % number of particles
wMax=1000; % max expected weight of a particle
w=rand(np,1)*wMax; % weight of particles

% add many light particles
np=1000; % number of particles
wMax=10; % max expected weight of a particle
w=[w;rand(np,1)*wMax]; % weight of particles


[ws,I]=sort(w,1,'descend');

for k=1:length(ws)
mMin=find(wM==min(wM),1);
wM(mMin)=wM(mMin)+ws(k);
figure(1)
bar([1:m],wM)
title([num2str(k) 'of ' num2str(length(ws))])
drawnow
end

%figure(1)
%bar([1:m],wM)