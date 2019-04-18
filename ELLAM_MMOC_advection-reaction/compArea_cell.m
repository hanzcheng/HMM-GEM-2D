function [AreaIn,AreaT]=compArea_cell(ncell,pointsToTrack, newPoints,cell_ep)

    AreaIn=zeros(ncell,ncell); %AreaIn(i,j) gives the area of the intersection of cell i and tracked cell j
    origPoints=cell(1,ncell);
    trackedPoly=cell(1,ncell);
    AreaT=zeros(1,ncell);
    P1=[];
    P2=[];
    P3=[];
    for i=1:ncell
        ptsInvolved = cell_ep(i,:);
        ptsInvolved(ptsInvolved == 0)=[];
        origPoints{i}=pointsToTrack(ptsInvolved,:);
        trackedPoly{i}=newPoints(ptsInvolved,:);
        AreaT(i)=compute_area([trackedPoly{i}(:,1) trackedPoly{i}(:,2)]);
    end
     
    for i=1:ncell
        P1.x=origPoints{i}(:,1);
        P1.y=origPoints{i}(:,2);
        for j=1:ncell
            P2.x=trackedPoly{j}(:,1);
            P2.y=trackedPoly{j}(:,2);
            P3=PolygonClip(P1,P2,1);
            numPol=length(P3);
            for p=1:numPol
            P3(p).x=[P3(p).x; P3(p).x(1)];
            P3(p).y=[P3(p).y; P3(p).y(1)];
            vertInt=[P3(p).x P3(p).y];
            AreaIn(i,j)=AreaIn(i,j)+compute_area(vertInt);
            end
%             if i==1 && j==1
%                 Pol1=[P1.x P1.y]
%                 Pol2=[P2.x P2.y]
%                 Pol3=[P3.x P3.y]
%             end
        end
    end
%     AreaIn=sparse(AreaIn);
end