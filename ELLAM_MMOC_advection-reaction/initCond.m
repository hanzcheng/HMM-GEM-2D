function uv=initCond(z)
global t 
global testcase
x=z(1); 
y=z(2);
uv=zeros(size(z,1),1);
if (testcase==0)
	uv = ((x-.25).^2+(y-.75).^2 < 1/64);	
elseif testcase == 1
    uv(x>=0.5)=1;
elseif testcase == 2
uv=t^2*(x .^ 2 + y .^ 2);

elseif testcase == 3
uv=(x-x.^2).*(y-y.^2)*(t+t^2);

elseif testcase ==4
uv=sin(pi*x) .* sin(pi*y) .* sin(t);

elseif testcase == 5 %the rectangular region (0.2,0.4) x (0.1,0.3)
    nZx1=find(x>=0.2);
    nZx2=find(x<=0.4);
    nZx=intersect(nZx1,nZx2);
    nZy1=find(y>=0.1);
    nZy2=find(y<=0.3);
    nZy=intersect(nZy1,nZy2);
    nZ=intersect(nZx,nZy);
    
uv(nZ)=1;

%     nZx1a=find(x>=0.2);
%     nZx2a=find(x<=0.4);
%     nZxa=intersect(nZx1a,nZx2a);
%     nZy1a=find(y>=0.6);
%     nZy2a=find(y<=0.7);
%     nZya=intersect(nZy1a,nZy2a);
%     nZa=intersect(nZxa,nZya);
%     
% uv(nZa)=1;

end

