function uv=initCond(z)
global testcase
x=z(:,1); 
y=z(:,2);
uv=zeros(size(z,1),1);
if (testcase==1) % square block (0.0625,0.0625)x(0.3125,0.3125)
	
    nZx1=find(x>=0.0625);
    nZx2=find(x<=0.3125);
    nZx=intersect(nZx1,nZx2);
    nZy1=find(y>=0.0625);
    nZy2=find(y<=0.3125);
    nZy=intersect(nZy1,nZy2);
    nZ=intersect(nZx,nZy);
    
uv(nZ)=1;
    
elseif testcase == 2 
 uv = exp(-10*((x-.1875) .^ 2 + (y-.1875) .^ 2));
uv(x<0) = 0;

elseif testcase == 3
    uv = ((x-.25).^2+(y-.75).^2 < 1/64);

elseif testcase == 4
    uv = exp(-10*((x-.25) .^ 2 + (y-.75) .^ 2));
end

