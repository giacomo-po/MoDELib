clc
clear all
close all

Theta=[0 1 0.5]; % array of time integtaion parameters. 0=forward Euler, 1=backward Euler, 0.5=Crank-Nicolson
nSteps=5;        % number of time integration steps
D=0.0005;        % diffusion coefficient

fontSize=16;
ThetaClr={'-b','--r','-.g'};
ThetaNames={'forward Euler','backward Euler','Crank-Nicolson'}
comp={'x_1','value'};
clrCol=2;

figure(1)
hold on

for t=1:length(Theta)
    theta=Theta(t);
    system('rm F/F_*');
    system(['./ibvp ' num2str(1) ' ' num2str(theta) ' ' num2str(nSteps) ' ' num2str(D)])
    
    for stepID=nSteps
        
        S=load(['F/F_' num2str(stepID) '.txt']);
        
        sMax=max(S(:,clrCol));
        sMin=min(S(:,clrCol));
        clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));
        
        for e=1:2:size(S,1)
            H(t)=plot(S(e:e+1,1),S(e:e+1,2),ThetaClr{t},'Linewidth',2);
        end
        grid on
        xlabel('X','FontSize',fontSize)
        axis([-1 1 -0.2 1.2])
        legend(H,ThetaNames)
        title(['D=' num2str(D) ', nSteps=' num2str(nSteps) ', dt=1/nSteps'])
        set(gca,'FontSize',fontSize)
        drawnow
    end
    
end
