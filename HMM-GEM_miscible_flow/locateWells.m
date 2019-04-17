function [icell,pcell,iCellCount,pCellCount]= locateWells(xprod,yprod,xinj,yinj,vertex,cell_v,ncell,thres)
nPts = length(xprod); %check how many injection/production points (should be equal)
pCellCount=0; %counts the number of cells where the production and injection points are located
iCellCount=0;
icell=[];
pcell=[];
for j=1:nPts
    ppoint=find(abs(vertex(:,1)-xprod(j))<thres & abs(vertex(:,2)-yprod(j))<thres);% production point
    ipoint=find(abs(vertex(:,1)-xinj(j))<thres & abs(vertex(:,2)-yinj(j))<thres);% injection point
    for i=1:ncell %loop over the cells to locate the injection and production cells
        if find(cell_v{i} == ppoint)
            pCellCount=pCellCount+1;
            pcell(pCellCount)=i;
        elseif find(cell_v{i} == ipoint)
            iCellCount=iCellCount+1;
            icell(iCellCount)=i;
        end
    end
end
end