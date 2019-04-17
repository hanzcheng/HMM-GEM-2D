function ec = edge_centers(currcell,cell_e,cell_v,vertex)
%% returns the centers of the edges of the current cell
ne=size(cell_e{currcell},2);
ec=zeros(ne,2);
for i=1:ne
    ec(i,:)= 0.5* (vertex(cell_v{currcell}(i),:)+ vertex(cell_v{currcell}(i+1),:));
end
end