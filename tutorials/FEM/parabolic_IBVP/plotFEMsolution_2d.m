clc
clear all
close all

Theta=[0 1 0.5]; % array of time integtaion parameters. 0=forward Euler, 1=backward Euler, 0.5=Crank-Nicolson
nSteps=3;        % number of time integration steps
D=0.001;        % diffusion coefficient

fontSize=16;
%ThetaClr={'-b','--r','-.g'};
ThetaNames={'forward Euler','backward Euler','Crank-Nicolson'}
comp={'x_1','x_2','value'};
clrCol=3;

for t=1:length(Theta)
    theta=Theta(t);
    
    system('rm F/F_*');
    system(['./ibvp ' num2str(2) ' ' num2str(theta) ' ' num2str(nSteps) ' ' num2str(D)])

for stepID=nSteps
    
        figure(t)
         hold on

        S=load(['F/F_' num2str(stepID) '.txt']);
        
        
        sMax=max(S(:,clrCol))
        sMin=min(S(:,clrCol))
        clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));
        for f=1:3:size(S,1)
            fill3(S(f:f+2,1),S(f:f+2,2),S(f:f+2,3),clrID(f:f+2));
        end
        grid on
        xlabel('X','FontSize',fontSize)
        ylabel('Y','FontSize',fontSize)
        legend(ThetaNames{t})
        title(['D=' num2str(D) ', nSteps=' num2str(nSteps) ', dt=1/nSteps'])
        set(gca,'FontSize',fontSize)
        set(gca,'View',[40 40])
    drawnow
end

end



