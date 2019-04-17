function sub_cell_ep=whichPoint(ncell,sub_ncell,cell_e,nxEdges,nyEdges,pointsToTrack,MpPerEdge)
%later on, add pPerEdge to the parameters needed
ctr=1;
sub_cell_ep=zeros(sub_ncell,MpPerEdge);
nbPts=size(pointsToTrack,1);
tol=1e-8;
for i=1:ncell
    nbe=length(cell_e{i});
    for j=1:nbe
        for s=1:size(nxEdges{i}{j},1)
            xydiff=[nxEdges{i}{j}(s)*ones(nbPts,1) nyEdges{i}{j}(s)*ones(nbPts,1)]-pointsToTrack;
            normXy= sqrt(sum(abs(xydiff).^2,2));
            ptLoc=find(normXy<=tol);
            sub_cell_ep(ctr,s)=ptLoc;
        end
        ctr=ctr+1;
    end
end