function [changedArea,pctErr,aErr,nAdj,pctErrTbInit] = adjust_volumes_for_local_mass(AreaInF,areaToChange,alffa,invTF,area,ncell,cell_e,cell_n,newPoints,nextTri,triInvI,subCellsIn,sub_cell_ep,mainCell,pOnEdge,icell,pcell,KRDarcyU,varargin)
if ~isempty(varargin)
    cellsInv = varargin{1};
    cellsUnInv = varargin{2}; %cells uninvolved with current process
    startGEM = varargin{3};
    pNear = varargin{4};
    isELLAM = varargin{5}; % 1 is ellam 0 is mmoc
else
    cellsInv = 1:ncell;
    startGEM = 0;
end

nCellsInv = length(cellsInv);
pctErr = zeros(ncell,1);
aErr = zeros(ncell,1);
eCells = [icell pcell]; % cells to exclude in the post-processing
nAdj=1;
iCellCount = length(icell);
changedArea = areaToChange;
for i=1:nCellsInv
    currCell = cellsInv(i);
    if ~any(eCells==currCell)
        pctErr(currCell) = abs(sum(changedArea(:,currCell))-area(currCell))/area(currCell);
        aErr(currCell) = sum(changedArea(:,currCell))-area(currCell); % if +, excess, if -, missing
    end
end
pctErrTbInit = pctErr;
aTF=zeros(iCellCount,1);
%% include regions near the  injection well in local mass cons (approximate only)
if startGEM==0
    for ji=1:iCellCount
        aTF(ji)=sum(AreaInF(invTF{ji},icell(ji)));
        
        for k=1:size(invTF{ji},1)
            aExpect = area(invTF{ji}(k))-AreaInF(invTF{ji}(k),icell(ji))+AreaInF(invTF{ji}(k),icell(ji))/(aTF(ji))...
                *(1-exp(-alffa))*area(invTF{ji}(k));
            pctErr(invTF{ji}(k)) = abs(sum(changedArea(:,invTF{ji}(k)))-aExpect)/aExpect;
            aErr(invTF{ji}(k)) = sum(changedArea(:,invTF{ji}(k)))-aExpect;
%%            uncomment these 2 lines if we dont want to include cells around injection for local mass (as in hexa)
%             pctErr(invTF{ji}(k)) = 0;
%             aErr(invTF{ji}(k)) = 0;
        end
    end
end
pctErr(eCells)=0;
while (max(pctErr)>1e-4 && nAdj<100)
    
    %% adjust areas to remove errors
    for i=1:nCellsInv
        currCell = cellsInv(i);
        if ~any(eCells==currCell)
            nbe = length(cell_e{currCell});
            eInv = zeros(nbe,4); % edge involved, and the cells to be modified (indexed at 2 for orig cell, 3 for tracked cell)
            % 4th index is for the area that will be affected (should
            % check that it wont become negative)
            vInv = zeros(nbe,2);
            vMagE = zeros(nbe,1); %magnitude of the velocity along edges
            vMagV = zeros(nbe,1); %magnitude of the velocity along vertices
            ePtCheck = zeros(nbe,2); %edge pts to check
            vPtCheck = zeros(nbe,2); %vertices to check
            for j=1:nbe % loop over the edges to detect where to add volume
                if cell_n{currCell}(j)>0 % only adjust edges which will affect "next" cell (if > i) otherwise adjust any, except boundary
                    nEdgePts = pOnEdge(cell_e{currCell}(j))+2; %pts along the edge (including the 2 vertices)
                    midEdge = ceil(nEdgePts/2);
                    locEdge = find(pNear==sub_cell_ep(subCellsIn{currCell}(j),midEdge));
                    ePtCheck(j,:) = newPoints(locEdge,:);
                    ptTri = nextTri(locEdge); % triangle for the point
                    velPt = KRDarcyU(1,ptTri)*ePtCheck(j,:)' + KRDarcyU(2:3,ptTri);%evaluate the velocity at the given pt
                    locV1 = find(pNear==sub_cell_ep(subCellsIn{currCell}(j),1));
                    Tv1 = newPoints(locV1,:); % tracked vertices
                    locV2 = find(pNear==sub_cell_ep(subCellsIn{currCell}(j),nEdgePts));
                    Tv2 = newPoints(locV2,:); % tracked vertices
                    oN = (Tv2 - Tv1)/norm(Tv2 - Tv1)*[0 1; -1 0]; % outward normal
                    uDotN = (1)*oN*velPt;
                    if startGEM==0
                        uDotN = (-1)^(nAdj+1)*uDotN;
                    end
                    if aErr(currCell)>0
                        if uDotN<0
                            eInv(j,1) = 1;
                            eInv(j,2) = mainCell(ptTri);
                            eInv(j,3) = cell_n{currCell}(j);
                            eInv(j,4) = changedArea(eInv(j,2),eInv(j,3));
                            vMagE(j) = norm(velPt);
                        end
                    elseif aErr(currCell)<0
                        if uDotN>0
                            eInv(j,1) = 1;
                            eInv(j,2) = mainCell(ptTri);
                            eInv(j,3) = cell_n{currCell}(j);
                            eInv(j,4) = changedArea(eInv(j,2),eInv(j,3));
                            vMagE(j) = norm(velPt);
                        end
                    end
                    if j>1
                        if (eInv(j-1,1)==1 && eInv(j,1)==1)
                            vInv(j,1) = 1;
                            vPtCheck(j,:) = Tv1;
                            vertTri = nextTri(locV1);%location of vertex
                            velVert = KRDarcyU(1,vertTri)*Tv1' + KRDarcyU(2:3,vertTri);%evaluate the velocity at the given pt
                            vMagV(j) = norm(velVert);
                            vInv(j,2) = mainCell(vertTri);
                        end
                    end
                end
            end
            
            if (eInv(nbe,1)==1 && eInv(1,1)==1)
                vInv(1,1) = 1;
                locV1 = find(pNear == sub_cell_ep(subCellsIn{currCell}(1),1));
                vPtCheck(1,:) = newPoints(locV1,:); % tracked vertex
                ptTri = nextTri(locV1);
                velPt = KRDarcyU(1,ptTri)*vPtCheck(1,:)' + KRDarcyU(2:3,ptTri);%evaluate the velocity at the given pt
                vMagV(1) = norm(velPt);
                vInv(1,2) = mainCell(ptTri);
            end
            totalMag = sum(vMagV) + sum(vMagE);
            %distribute a proper percentage of the area to the cells
            %around F_-dt(K). Also detect which F_-dt(Ki) are affected
            for j=1:nbe
                vInv(j,1) = vMagV(j)/totalMag;
                eInv(j,1) = vMagE(j)/totalMag;
                if eInv(j,2)~=0
                    if changedArea(eInv(j,2),currCell) - eInv(j,1)*aErr(currCell)>0
                        changedArea(eInv(j,2),currCell) = changedArea(eInv(j,2),currCell) - eInv(j,1)*aErr(currCell);
                        if isELLAM
                            finTest = eInv(j,2);
                        else
                            finTest = eInv(j,3);
                        end
                        if (changedArea(eInv(j,2),eInv(j,3)) +  eInv(j,1)*aErr(currCell)>0 || any(cellsUnInv==finTest))
                            changedArea(eInv(j,2),eInv(j,3)) =  changedArea(eInv(j,2),eInv(j,3)) +  eInv(j,1)*aErr(currCell);
                            aErr(eInv(j,3)) = aErr(eInv(j,3))+ eInv(j,1)*aErr(currCell);
                        else
                            aErr(eInv(j,3)) = aErr(eInv(j,3))- changedArea(eInv(j,2),eInv(j,3));
                            changedArea(eInv(j,2),currCell) = changedArea(eInv(j,2),currCell) + eInv(j,1)*aErr(currCell) + changedArea(eInv(j,2),eInv(j,3));
                            changedArea(eInv(j,2),eInv(j,3)) =  0;
                        end
                    else
                            changedArea(eInv(j,2),eInv(j,3)) =  changedArea(eInv(j,2),eInv(j,3)) +  changedArea(eInv(j,2),currCell);
                            aErr(eInv(j,3)) = aErr(eInv(j,3))+ changedArea(eInv(j,2),currCell);
                            changedArea(eInv(j,2),currCell) = 0;
                    end
                    
                end
                if vInv(j,2)~=0
                    locV1 = find(pNear == sub_cell_ep(subCellsIn{currCell}(j),1));
                    vertCells = mainCell(triInvI(locV1,triInvI(locV1,:)~=0));
                    vertCells(vertCells == currCell)=[];
                    vertCells = unique(vertCells);
                    nVertCells = length(vertCells);
                    if changedArea(vInv(j,2),currCell) - vInv(j,1)*aErr(currCell)>0
                        changedArea(vInv(j,2),currCell) = changedArea(vInv(j,2),currCell) - vInv(j,1)*aErr(currCell);
                        for v=1:nVertCells
                        if isELLAM
                            finTest = vInv(j,2);
                        else
                            finTest = vertCells(v);
                        end
                            if (changedArea(vInv(j,2),vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells>0 || any(cellsUnInv==finTest))
                                changedArea(vInv(j,2),vertCells(v)) = changedArea(vInv(j,2),vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells;
                                aErr(vertCells(v)) = aErr(vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells;
                            else
                                changedArea(vInv(j,2),currCell) = changedArea(vInv(j,2),currCell) + vInv(j,1)*aErr(currCell)/nVertCells+changedArea(vInv(j,2),vertCells(v));
                                aErr(vertCells(v)) = aErr(vertCells(v))-changedArea(vInv(j,2),vertCells(v));
                                changedArea(vInv(j,2),vertCells(v)) = 0;
                            end
                        end
                    else
                                changedArea(vInv(j,2),vertCells) = changedArea(vInv(j,2),vertCells)+changedArea(vInv(j,2),currCell)/nVertCells;
                                aErr(vertCells) = aErr(vertCells)+changedArea(vInv(j,2),currCell)/nVertCells;
                                changedArea(vInv(j,2),currCell) = 0;
                    end
                end
                
            end
            % subtract the amt added to the respective cells involved
        end
    end
    nAdj=nAdj+1;
    for i=1:nCellsInv
        currCell = cellsInv(i);
        if ~any(eCells==currCell)
            pctErr(currCell) = abs(sum(changedArea(:,currCell))-area(currCell))/area(currCell);
            aErr(currCell) = sum(changedArea(:,currCell))-area(currCell); % if +, excess, if -, missing
        end
    end
    if startGEM==0
        for ji=1:iCellCount
            for k=1:size(invTF{ji},1)
                aExpect = area(invTF{ji}(k))-AreaInF(invTF{ji}(k),icell(ji))+AreaInF(invTF{ji}(k),icell(ji))/(aTF(ji))...
                    *(1-exp(-alffa))*area(invTF{ji}(k));
                pctErr(invTF{ji}(k)) = abs(sum(changedArea(:,invTF{ji}(k)))-aExpect)/aExpect;
                aErr(invTF{ji}(k)) = sum(changedArea(:,invTF{ji}(k)))-aExpect;
                %% uncomment these 2 lines if we dont want to include cells around injection for local mass (as in hexa)
%                             pctErr(invTF{ji}(k)) = 0;
%             aErr(invTF{ji}(k)) = 0;
            end
        end
    end
    
    pctErr(eCells)=0;
end
end