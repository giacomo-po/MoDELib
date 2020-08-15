clear all
%close all
clc

frameID=5;
d=load(['H/H_' num2str(frameID) '.txt']);

%x=d(:,3)./(d(:,1)+d(:,2));
x=d(:,3)./(d(:,1)+d(:,2));
xs=sort(x);
ids=find(diff(xs));

x1=d(:,3)./(d(:,1));
x1s=sort(x1);
ids1=find(diff(x1s));

%x2=d(:,3)./(d(:,2));
%x2s=sort(x2);
%ids2=find(diff(x2s));


% Xu=unique(x);
% Lu=zeros(length(Xu),1);
% for i=1:length(Xu)
% Lu(i)=length(find(x<=Xu(i)));
% end

%%

figure
clf
hold on
plot(xs(ids),ids/length(xs),'Linewidth',2)
plot(x1s(ids1),ids1/length(x1s),'--','Linewidth',2)
%set(gca,'XScale','log')
grid on

x2=d(:,1);
x2s=sort(x2);
ids2=find(diff(x2s));

figure
clf
hold on
plot(x2s(ids2),ids2/length(x2s),'Linewidth',2)
%plot(x2s(ids2),ids2/length(x2s),'--','Linewidth',2)
%set(gca,'XScale','log')
grid on

