close all
clear all
clc

fontSize=16;
SIGMA=[3 5]

for s=1:length(SIGMA)
    sigma=SIGMA(s);
    
    %fileName=['SimpleCubic_sigma_' num2str(sigma)];
    fileName=['FCC_sigma_' num2str(sigma)];
    %fileName=['BCC_sigma_' num2str(sigma)];

    M=load([fileName '.txt']);
    A=M(1:3,:);
    B=M(4:6,:);
    CSL=M(7:9,:);
    DSCL=M(10:12,:);
    
    L=3.3;
    figure(2*s)
    clf
    
    plotLattice(DSCL,L,'k.',1)
    plotLattice(CSL,L,'g.',2)    
    plotLattice(A,L,'bo',1)
    plotLattice(B,L,'rx',1)
    legend('DSCL','CSL','A','B')
    title(strrep(fileName,'_',' '),'FontSize',fontSize)
    set(gca,'FontSize',fontSize)
    axis equal
    
end

