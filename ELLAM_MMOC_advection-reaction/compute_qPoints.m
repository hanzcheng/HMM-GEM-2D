
function [cB,mainCell] = compute_qPoints(ncell,cell_v,vertex,snBalls)
nBalls=snBalls^2;
bTotal=ncell*nBalls;
cB=zeros(bTotal,2); %ball center
mainCell=zeros(bTotal,1);
ctr=1;
for i=1:ncell
    xL=min(vertex(cell_v{i},1));
    xR=max(vertex(cell_v{i},1));
    yL=min(vertex(cell_v{i},2));
    yR=max(vertex(cell_v{i},2));
   
    xNow=linspace(xL+(xR-xL)/(2*snBalls),xR-(xR-xL)/(2*snBalls),snBalls);
    yNow=linspace(yL+(yR-yL)/(2*snBalls),yR-(yR-yL)/(2*snBalls),snBalls);
    [x,y]=meshgrid(xNow,yNow);
    x=x(:);
    y=y(:);
    cB(ctr:ctr+nBalls-1,1)=x;
    cB(ctr:ctr+nBalls-1,2)=y;
    mainCell(ctr:ctr+nBalls-1)=i;

    ctr=ctr+nBalls;
end


