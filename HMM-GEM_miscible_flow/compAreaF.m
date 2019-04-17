function [AreaIn,AreaT]=compAreaF(ncell,icell,subCellsIn,nxEdges,pointsToTrack, newPoints,sub_cell_ep,pFromInj)
    AreaIn=zeros(ncell,ncell); %AreaIn(i,j) gives the area of the intersection of cell i and tracked cell j
    origPoints=cell(1,ncell);
    tnewPoints=zeros(size(pointsToTrack));
    tnewPoints(pFromInj,:)=newPoints;
    newPoints=tnewPoints;
    trackedPoly=cell(1,length(icell));
    AreaT=zeros(1,ncell);
        sizeCell=zeros(1,ncell);
    for i=1:ncell
    nbe=size(subCellsIn{i},2);
    for j=1:nbe
        sizeCell(i)=sizeCell(i)+size(nxEdges{i}{j},1);
    end
    end
    for i=1:ncell
        nbe=size(subCellsIn{i},2);
          origPoints{i}=zeros(sizeCell(i),2);
        trackedPoly{i}=zeros(sizeCell(i),2);
        ctr=1;
        for j=1:nbe
             currSubcell=sub_cell_ep(subCellsIn{i}(j),:);
            currSubcell(currSubcell==0)=[];
            if any(icell==i)
            trackedPoly{i}(ctr:ctr+length(currSubcell)-1,:)=newPoints(currSubcell,:);
            end
            origPoints{i}(ctr:ctr+length(currSubcell)-1,:)=pointsToTrack(currSubcell,:);
            ctr=ctr+length(currSubcell);
        end        
        AreaT(i)=compute_area([trackedPoly{i}(:,1) trackedPoly{i}(:,2)]);
    end
    for i=1:ncell
        P1.x=origPoints{i}(:,1);
        P1.y=origPoints{i}(:,2);
        for j=1:length(icell)
            P2.x=trackedPoly{icell(j)}(:,1);
            P2.y=trackedPoly{icell(j)}(:,2);
            P3=PolygonClip(P1,P2,1);
            numPol=length(P3);
            for p=1:numPol
            P3(p).x=[P3(p).x; P3(p).x(1)];
            P3(p).y=[P3(p).y; P3(p).y(1)];
            vertInt=[P3(p).x P3(p).y];
            AreaIn(i,icell(j))=AreaIn(i,icell(j))+compute_area(vertInt);
            end
        end
    end
%     AreaIn=sparse(AreaIn);
end