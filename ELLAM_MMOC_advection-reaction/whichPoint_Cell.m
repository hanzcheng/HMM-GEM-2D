function cell_ep=whichPoint_Cell(ncell,cell_e,nxEdges,nyEdges,pointsToTrack,MpPerEdge)

cell_ep=zeros(ncell,MpPerEdge*4);
nbPts=size(pointsToTrack,1);
tol=1e-8;
for i=1:ncell
    ctr=1;
    nbe=length(cell_e{i});
    for j=1:nbe
        for s=1:size(nxEdges{i}{j},1)
            xydiff=[nxEdges{i}{j}(s)*ones(nbPts,1) nyEdges{i}{j}(s)*ones(nbPts,1)]-pointsToTrack;
            normXy= sqrt(sum(abs(xydiff).^2,2));
            ptLoc=find(normXy<=tol);
            
            cell_ep(i,ctr+s-1)=ptLoc;
        end
        ctr=ctr+size(nxEdges{i}{j},1);
    end
end