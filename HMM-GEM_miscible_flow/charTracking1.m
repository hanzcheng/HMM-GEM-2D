function [tpoints]=charTracking1(pointsToTrack,dtTrack,currTri,KRDarcyU,track)
%% Performs the Characteristic Tracking for RT0 velocity fields
x=pointsToTrack(:,1);
y=pointsToTrack(:,2);
xnew=zeros(length(x),1);
ynew=zeros(length(y),1);
a=KRDarcyU(1,currTri);
b1=KRDarcyU(2,currTri);
b2=KRDarcyU(3,currTri);
%here, KRDarcyU is of the form ax+[b1;b2]
a=a';
b1=b1';
b2=b2';
if track==0 %backward tracking
    xnew(a==0)=-b1(a==0) .* dtTrack(a==0)+x(a==0);
    ynew(a==0)=-b2(a==0) .* dtTrack(a==0)+y(a==0);
    xnew(a~=0)=((a(a~=0) .* x(a~=0) + b1(a~=0)) .* exp(-a(a~=0) .* dtTrack(a~=0))- b1(a~=0)) ./ a(a~=0);
    ynew(a~=0)=((a(a~=0) .* y(a~=0) + b2(a~=0)) .* exp(-a(a~=0) .* dtTrack(a~=0))- b2(a~=0)) ./ a(a~=0);
else        %forward tracking
    xnew(a==0)=b1(a==0) .* dtTrack(a==0)+x(a==0);
    ynew(a==0)=b2(a==0) .* dtTrack(a==0)+y(a==0);
    xnew(a~=0)=((a(a~=0) .* x(a~=0) + b1(a~=0)) .* exp(a(a~=0) .* dtTrack(a~=0))- b1(a~=0)) ./ a(a~=0);
    ynew(a~=0)=((a(a~=0) .* y(a~=0) + b2(a~=0)) .* exp(a(a~=0) .* dtTrack(a~=0))- b2(a~=0)) ./ a(a~=0);
end

%% fixes so that computational/machine error doesn't lead to something outside the domain
if pointsToTrack(1,1)==0 && pointsToTrack(1,2)==0
    xnew(1)=0;
    ynew(1)=0;
end
xO=find(xnew>1000);
xO1=find(xnew<1001);
xO=intersect(xO,xO1); 
yO=find(ynew>1000);
yO1=find(ynew<1001);
yO=intersect(yO,yO1);
xnew(xO)=1000;
ynew(yO)=1000;
tpoints=[xnew, ynew];
