function [xEdges,yEdges]= edge_points(ncell,cell_e,cell_v,vertex,nPts)
%% returns the x and y coordinates of nPts equally spaced points along the edges of the cells
% together with the edge vertices
%Note: if nPts=1 then we are simply looking at the edge centers
%To get the x coordinates of the points on cell i along the jth edge, we
%use xEdges{i}(:,j);
%% Note: this needs to be modified in Peaceman case to accumulate more traceback points if near injection
xEdges=cell(1,ncell);
yEdges=cell(1,ncell);
for i=1:ncell
    mult=0:nPts(i)+1;
    mult=mult';
    nbe=length(cell_e{i});
    mx=vertex(cell_v{i}(2:nbe+1),1)-vertex(cell_v{i}(1:nbe),1);
    my=vertex(cell_v{i}(2:nbe+1),2)-vertex(cell_v{i}(1:nbe),2);
    xe=repmat(vertex(cell_v{i}(1:nbe),1)',nPts(i)+2,1)+mult/(nPts(i)+1)*mx';
    ye=repmat(vertex(cell_v{i}(1:nbe),2)',nPts(i)+2,1)+mult/(nPts(i)+1)*my';
xEdges{i}=xe;
yEdges{i}=ye;
end
