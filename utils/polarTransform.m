function polartrans = polarTransform2( t , rstep, thetastep_deg )
%polarTransform Summary of this function goes here
%   Detailed explanation goes here
[N,~] = size(t);
t = double(t);

center = N/2 + 1;

thetastep = thetastep_deg*pi/180;

polartrans = zeros(N/2/rstep,360/thetastep_deg);

rs = 0:rstep:(N/2-rstep);
qs = 0:thetastep:(2*pi-thetastep);

[Rs,Qs] = meshgrid(rs,qs);

Ylocs = center + Rs.*sin(Qs);
Xlocs = center + Rs.*cos(Qs);

Xlocs = round(Xlocs);
Ylocs = round(Ylocs);
Xlocs(Xlocs == N + 1) = N;
Ylocs(Ylocs == N + 1) = N;

for i = 1:numel(rs)
    for j = 1:numel(qs)
        polartrans(i,j) = t(Ylocs(j,i),Xlocs(j,i));
    end
end



end

