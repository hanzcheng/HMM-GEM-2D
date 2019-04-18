% Computes the real gravity center of the cells, by splitting them into
% triangles
function cg = gravity_centers(ncell,cell_v,vertex,area);

cg=zeros(ncell,2);

for i=1:ncell
	vertex_loc=vertex(cell_v{i},:);
	% Nb of edges in the cell
	nbe = size(vertex_loc,1)-1;

	% We compute one point inside the cell, to split the cells into triangles
	ptK=sum(vertex_loc([1:nbe],:))/nbe;
			
	% Loop over vertices to create the triangles (ptK,vertex1,vertex2)
	for j=1:nbe
		% vertices of the triangle
		v1=vertex_loc(j,:);
		v2=vertex_loc(j+1,:);
		% area of triangle
		areatri=0.5*det([v1-ptK;v2-ptK]);
		% add the center of (ptK,v1,v2) ponderated by areatri
		cg(i,:)=cg(i,:)+areatri*(ptK+v1+v2)/3;
	end
	% Apply total weight
	cg(i,:)=cg(i,:)/area(i);
end;


