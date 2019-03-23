function c=scalar2color(v,vmin,vmax)

c=ones(1,3); % white
if (v < vmin)
    v = vmin;
end
if (v > vmax)
    v = vmax;
end
dv = vmax - vmin;


% c1=[0 0 1];
% c2=[0 1 1];
% c3=[0 1 0];
% c4=[1 1 0];
% c5=[1 0 0];

 c1=rgb('RoyalBlue');
 c2=rgb('LightBlue');
 c3=rgb('LightGreen');
 c4=rgb('Khaki');
 c5=rgb('Salmon');

% c5=rgb('RoyalBlue');
% c4=rgb('LightBlue');
% c3=rgb('LightGreen');
% c2=rgb('Khaki');
% c1=rgb('Salmon');

if (v < (vmin + 0.25 * dv))
    vLow =vmin;
    vHigh=vmin+0.25 * dv;
    cLow=c1;
    cHigh=c2;
elseif (v < (vmin + 0.5 * dv))
    vLow =vmin+0.25 * dv;
    vHigh=vmin+0.50 * dv;
    cLow=c2;
    cHigh=c3;
elseif (v < (vmin + 0.75 * dv))
    vLow =vmin+0.50 * dv;
    vHigh=vmin+0.75 * dv;
    cLow=c3;
    cHigh=c4;
else
    vLow =vmin+0.75 * dv;
    vHigh=vmin+1.00 * dv;
    cLow=c4;
    cHigh=c5;
end

c=(v-vLow)/(vHigh-vLow)*cHigh+(vHigh-v)/(vHigh-vLow)*cLow;

%    if (v < (vmin + 0.25 * dv))
%       c(1) = 0;
%       c(2) = 4 * (v - vmin) / dv;
%    elseif (v < (vmin + 0.5 * dv))
%       c(1) = 0;
%       c(3) = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
%    elseif (v < (vmin + 0.75 * dv))
%       c(1) = 4 * (v - vmin - 0.5 * dv) / dv;
%       c(3) = 0;
%    else
%       c(2) = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
%       c(3) = 0;
%    end

