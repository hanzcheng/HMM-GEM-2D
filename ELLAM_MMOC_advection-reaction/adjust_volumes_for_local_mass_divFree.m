function [changedArea,pctErr,aErr,nAdj,aErrTbInit] = adjust_volumes_for_local_mass_divFree(vertex,cell_v,areaToChange,area,ncell,cell_e,cell_n,newPointsB,cellInvI,cell_ep,pOnEdge,dv_exact)
cellsInv = 1:ncell;
nCellsInv = length(cellsInv);
pctErr = zeros(ncell,1);
aErr = zeros(ncell,1);
nAdj=1;

changedArea = areaToChange;
for i=1:nCellsInv
    currCell = cellsInv(i);
    pctErr(currCell) = abs(sum(changedArea(:,currCell))-area(currCell))/area(currCell);
    aErr(currCell) = sum(changedArea(:,currCell))-area(currCell); % if +, excess, if -, missing
end
aErrTbInit = aErr;
while (max(pctErr)>1e-6 && nAdj<ncell)
    
    %% adjust areas to remove errors
    for i=1:nCellsInv
        currCell = cellsInv(i);
        cellsInt = find(changedArea(:,currCell)>0); %cells intersected
        nbe = length(cell_e{currCell});
        eInv = zeros(nbe,4); % edge involved, and the cells to be modified (indexed at 2 for orig cell, 3 for tracked cell)
        % 4th index is for the area that will be affected (should
        % check that it wont become negative)
        vInv = zeros(nbe,2);
        vMagE = zeros(nbe,1); %magnitude of the velocity along edges
        vMagV = zeros(nbe,1); %magnitude of the velocity along vertices
        ePtCheck = zeros(nbe,2); %edge pts to check
        vPtCheck = zeros(nbe,2); %vertices to check
        eCtr = 0;
        for j=1:nbe % loop over the edges to detect where to add volume
            nEdgePts = pOnEdge(cell_e{currCell}(j))+2; %pts along the edge (including the 2 vertices)
            if cell_n{currCell}(j)>0 % only adjust edges which will affect "next" cell (if > i) otherwise adjust any, except boundary
                midEdge = ceil(nEdgePts/2);
                ePtCheck(j,:) = newPointsB(cell_ep(currCell,midEdge+eCtr),:);
                for s=1:length(cellsInt)
                    vertexNow = vertex(cell_v{cellsInt(s)},:);
                    xP = vertexNow(:,1);
                    yP = vertexNow(:,2);
                    testPoly = inpolygon(ePtCheck(j,1),ePtCheck(j,2),xP,yP);
                    if testPoly==1
                        ptCell = cellsInt(s);
                    end
                end
                velPt = dv_exact(ePtCheck(j,1),ePtCheck(j,2));
                Tv1 = newPointsB(cell_ep(currCell,1+eCtr),:); % tracked vertices
                Tv2 = newPointsB(cell_ep(currCell,nEdgePts+eCtr),:); % tracked vertices
                oN = (Tv2 - Tv1)/norm(Tv2 - Tv1)*[0 1; -1 0]; % outward normal
                uDotN = (-1)^(nAdj+1)*oN*velPt';
                if aErr(currCell)>0
                    if uDotN<0
                        eInv(j,1) = 1;
                        eInv(j,2) = ptCell;
                        eInv(j,3) = cell_n{currCell}(j);
                        eInv(j,4) = changedArea(eInv(j,2),eInv(j,3));
                        vMagE(j) = norm(velPt);
                    end
                elseif aErr(currCell)<0
                    if uDotN>0
                        eInv(j,1) = 1;
                        eInv(j,2) = ptCell;
                        eInv(j,3) = cell_n{currCell}(j);
                        eInv(j,4) = changedArea(eInv(j,2),eInv(j,3));
                        vMagE(j) = norm(velPt);
                    end
                end
                if j>1
                    if (eInv(j-1,1)==1 && eInv(j,1)==1)
                        vInv(j,1) = 1;
                        vPtCheck(j,:) = Tv1;
                        for s=1:length(cellsInt)
                            vertexNow = vertex(cell_v{cellsInt(s)},:);
                            xP = vertexNow(:,1);
                            yP = vertexNow(:,2);
                            testPoly = inpolygon(Tv1(1),Tv1(2),xP,yP);
                            if testPoly==1
                                ptCell = cellsInt(s);
                            end
                        end
                        velVert = dv_exact(Tv1(1),Tv1(2)); %evaluate the velocity at the given pt
                        vMagV(j) = norm(velVert);
                        vInv(j,2) = ptCell;
                    end
                end
            end
            eCtr = eCtr+nEdgePts;
        end
        
        if (eInv(nbe,1)==1 && eInv(1,1)==1)
            vInv(1,1) = 1;
            vPtCheck(1,:) = newPointsB(cell_ep(currCell,1),:); % tracked vertex
            velPt = dv_exact(vPtCheck(1,1),vPtCheck(1,2)); %evaluate the velocity at the given pt
            vMagV(1) = norm(velPt);
            for s=1:length(cellsInt)
                vertexNow = vertex(cell_v{cellsInt(s)},:);
                xP = vertexNow(:,1);
                yP = vertexNow(:,2);
                testPoly = inpolygon(vPtCheck(j,1),vPtCheck(j,2),xP,yP);
                if testPoly==1
                    ptCell = cellsInt(s);
                end
            end
            vInv(1,2) = ptCell;
        end
        totalMag = sum(vMagV) + sum(vMagE);
        %distribute a proper percentage of the area to the cells
        %around F_-dt(K). Also detect which F_-dt(Ki) are affected
        eCtr = 0;
        for j=1:nbe
            nEdgePts = pOnEdge(cell_e{currCell}(j))+2; %pts along the edge (including the 2 vertices)
            vInv(j,1) = vMagV(j)/totalMag;
            eInv(j,1) = vMagE(j)/totalMag;
            if eInv(j,2)~=0
                if changedArea(eInv(j,2),currCell) - eInv(j,1)*aErr(currCell)>0
                    changedArea(eInv(j,2),currCell) = changedArea(eInv(j,2),currCell) - eInv(j,1)*aErr(currCell);
                    if changedArea(eInv(j,2),eInv(j,3)) +  eInv(j,1)*aErr(currCell)>0 % after adjustment, cell should have + volume
                        changedArea(eInv(j,2),eInv(j,3)) =  changedArea(eInv(j,2),eInv(j,3)) +  eInv(j,1)*aErr(currCell);
                        aErr(eInv(j,3)) = aErr(eInv(j,3))+ eInv(j,1)*aErr(currCell);
                    else
                        aErr(eInv(j,3)) = aErr(eInv(j,3)) - changedArea(eInv(j,2),eInv(j,3));
                        changedArea(eInv(j,2),currCell) = changedArea(eInv(j,2),currCell) + eInv(j,1)*aErr(currCell) + changedArea(eInv(j,2),eInv(j,3));
                        changedArea(eInv(j,2),eInv(j,3)) = 0;
                    end
                else
                    if changedArea(eInv(j,2),eInv(j,3)) +  changedArea(eInv(j,2),currCell) < area(eInv(j,2))
                        changedArea(eInv(j,2),eInv(j,3)) =  changedArea(eInv(j,2),eInv(j,3)) +  changedArea(eInv(j,2),currCell);
                        aErr(eInv(j,3)) = aErr(eInv(j,3))+ changedArea(eInv(j,2),currCell);
                        changedArea(eInv(j,2),currCell) = 0;
                    else
                        %                                 aErrTb(eInv(j,3)) = aErrTb(eInv(j,3)) - AreaInB(eInv(j,2),eInv(j,3))+area(eInv(j,2));
                        %                                 AreaInB(eInv(j,2),i) = AreaInB(eInv(j,2),i) + AreaInB(eInv(j,2),eInv(j,3))-area(eInv(j,2));
                        %                                 AreaInB(eInv(j,2),eInv(j,3)) = area(eInv(j,2));
                    end
                end
            end
            if vInv(j,2)~=0
                vertCells = cellInvI(cell_ep(currCell,1+eCtr),cellInvI(cell_ep(currCell,1+eCtr),:)~=0);
                vertCells(vertCells == currCell)=[];
                vertCells = unique(vertCells);
                nVertCells = length(vertCells);
                if changedArea(vInv(j,2),currCell) - vInv(j,1)*aErr(currCell)>0
                    changedArea(vInv(j,2),currCell) = changedArea(vInv(j,2),currCell) - vInv(j,1)*aErr(currCell);
                    for v=1:nVertCells
                        if changedArea(vInv(j,2),vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells>0
                            changedArea(vInv(j,2),vertCells(v)) = changedArea(vInv(j,2),vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells;
                            aErr(vertCells(v)) = aErr(vertCells(v))+vInv(j,1)*aErr(currCell)/nVertCells;
                        else
                            aErr(vertCells(v)) = aErr(vertCells(v))-changedArea(vInv(j,2),vertCells(v));
                            changedArea(vInv(j,2),currCell) = changedArea(vInv(j,2),currCell) + vInv(j,1)*aErr(currCell)/nVertCells+changedArea(vInv(j,2),vertCells(v));
                            changedArea(vInv(j,2),vertCells(v)) = 0;
                        end
                    end
                else
                    for v=1:nVertCells
                        if changedArea(vInv(j,2),vertCells(v))+changedArea(vInv(j,2),currCell)/nVertCells<area(vInv(j,2))
                            changedArea(vInv(j,2),vertCells(v)) = changedArea(vInv(j,2),vertCells(v))+changedArea(vInv(j,2),currCell)/nVertCells;
                            aErr(vertCells(v)) = aErr(vertCells(v))+changedArea(vInv(j,2),currCell)/nVertCells;
                            changedArea(vInv(j,2),currCell) = changedArea(vInv(j,2),currCell) - changedArea(vInv(j,2),currCell)/nVertCells;
                        else
                            %                                     aErrTb(vertCells(v)) = aErrTb(vertCells(v))-AreaInB(vInv(j,2),vertCells(v))+area(vInv(j,2));
                            %                                     AreaInB(vInv(j,2),i) = AreaInB(vInv(j,2),i)+AreaInB(vInv(j,2),vertCells(v))-area(vInv(j,2));
                            %                                     AreaInB(vInv(j,2),vertCells(v)) = area(vInv(j,2));
                        end
                    end
                end
            end
            eCtr = eCtr+nEdgePts;
        end
        % subtract the amt added to the respective cells involved
    end
    nAdj=nAdj+1;
    for i=1:nCellsInv
        currCell = cellsInv(i);
        pctErr(currCell) = abs(sum(changedArea(:,currCell))-area(currCell))/area(currCell);
        aErr(currCell) = sum(changedArea(:,currCell))-area(currCell); % if +, excess, if -, missing
    end

end
end