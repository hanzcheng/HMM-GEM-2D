function [nearInjCell,nearProdCell,nearProdCellOnly,TriNearInj,TriNearProd,pNearInj,pNearProd]=classifyPtsandRegions(ncell,center,diam,mnbe,cell_e,sub_cell_ep,subCellsIn,icell,pcell)
allCells = 1:ncell;
currDist1=1000*sqrt(2)*ones(1,ncell); % distance to injection well
currDist2=1000*sqrt(2)*ones(1,ncell); % distance to production well
nearInjCell=zeros(1,ncell);
iCellCount = length(icell);
pCellCount = length(pcell);
for i=1:ncell
    for j=1:iCellCount
        currDist1(i)=min(currDist1(i),sqrt(sum(center(i,:)-center(icell(j),:)).^2));
    end
end
for i=1:ncell
    for j=1:pCellCount
        currDist2(i)=min(currDist2(i),sqrt(sum(center(i,:)-center(pcell(j),:)).^2));
    end
    
    if currDist1(i)<currDist2(i) %nearer to the injection well, trace forward
        nearInjCell(i)=i;
    end
end

nearProdCellOnly=allCells-nearInjCell;
nearProdCellOnly(nearProdCellOnly==0)=[];
diffDist = abs(currDist1 - currDist2);
nearBoth = find(diffDist<mean(diam));
nearProdCell = [nearProdCellOnly nearBoth];
nearProdCell(nearProdCell==0)=[];
nearProdCell = unique(nearProdCell);
nearInjCell(nearInjCell==0)=[];
%%
TriNearInj=zeros(1,mnbe*length(nearInjCell));%triangles near the injection cells
ctr=1;
for i=1:length(nearInjCell)
    nbe=length(cell_e{nearInjCell(i)});
    TriNearInj(ctr:ctr+nbe-1)=subCellsIn{nearInjCell(i)}; %to get the triangles involved with the injection well(s)
    ctr=ctr+nbe;
end
TriNearInj(TriNearInj==0)=[];
pNearInj=sub_cell_ep(TriNearInj,:);
pNearInj(pNearInj==0)=[];
pNearInj=unique(pNearInj(:));

%%
TriNearProd=zeros(1,mnbe*length(nearProdCell));%triangles near the production cells
ctr=1;
for i=1:length(nearProdCell)
    nbe=length(cell_e{nearProdCell(i)});
    TriNearProd(ctr:ctr+nbe-1)=subCellsIn{nearProdCell(i)}; %to get the triangles involved with the production well(s)
    ctr=ctr+nbe;
end
TriNearProd(TriNearProd==0)=[];
pNearProd=sub_cell_ep(TriNearProd,:);
pNearProd(pNearProd==0)=[];
pNearProd=unique(pNearProd(:));
end