% Write an ensight version of the solution
%
% uu should be a vector of size nvert
%
function out=write_solution_ensight(todo,uu,iter,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg)

% Have to split the cells into triangles
nsubtri = 0;
for i=1:ncell
	nsubtri = nsubtri + size(cell_e{i},2);
end


if (todo==0) 

	cd ensight;
	delete ensight.*;
	delete Un/*;

	% ensight.case file
	fid=fopen('ensight.case','w');

	fprintf(fid,'FORMAT\ntype: ensight gold\n\n');
	fprintf(fid,'GEOMETRY\nmodel: 1 ensight.geo\n\n');
	fprintf(fid,'VARIABLE\nscalar per node: 1 Un Un/Un*****\n\n');
	fprintf(fid,'TIME\ntime set: 1\n');
	fprintf(fid,'number of steps: %d\n',size(times,2));
	fprintf(fid,'filename start number: 0\nfilename increment: 1\n');
	fprintf(fid,'time values:\n');
	for i=1:size(times,2);
		fprintf(fid,'%E\n',times(i));
	end;

	fclose(fid);

	% ensight.geo
	fid=fopen('ensight.geo','w');

	fprintf(fid,'Code EF P1 Clement\nMaillage format ensight\nnode id assign\nelement id assign\npart\n 1\n');
	fprintf(fid,'AllNodes\ncoordinates\n%d\n',ncell+nvert);
	for i=1:ncell
		fprintf(fid,'%E\n',cg(i,1));
	end
	for i=1:nvert
		fprintf(fid,'%E\n',vertex(i,1));
	end;
	for i=1:ncell
		fprintf(fid,'%E\n',cg(i,2));
	end
	for i=1:nvert
		fprintf(fid,'%E\n',vertex(i,2));
	end;
	for i=1:ncell
		fprintf(fid,'%E\n',0);
	end;
	for i=1:nvert
		fprintf(fid,'%E\n',0);
	end;
	fprintf(fid,'tria3\n%d\n',nsubtri);
	for i=1:ncell
		for j=1:size(cell_e{i},2)
			fprintf(fid,'%d %d %d\n',i,ncell+cell_v{i}(j),ncell+cell_v{i}(j+1));
		end
	end

	fclose(fid);
	
	cd ../;

else


	% We first compute the value of the solution at the vertices.
	% The value at a vertex s is the average of the values in the cells K around s,
	% ponderated by a weitgh equal to the angle of the cell K at s.


	uvert=zeros(nvert,1);
	sumcoefvert=zeros(nvert,1);

	for i=1:ncell
		for j=1:size(cell_v{i},2)-1
			% The vertices around cell_v{i}(j)=j1 are j0 and j2 computed as follows
			j1=cell_v{i}(j);
			j2=cell_v{i}(j+1);
			if (j==1)
				j0=cell_v{i}(size(cell_v{i},2)-1);
			else
				j0=cell_v{i}(j-1);
			end
			% The angle of the cell i around j1 is arccos of the normalised scalar product
			% between j0j1 and j1j2
			vec2=vertex(j2,:)-vertex(j1,:);
			vec1=vertex(j0,:)-vertex(j1,:);
			scal=dot(vec1,vec2);
			% Cross-product between the vectors
			crossprod=vec2(1)*vec1(2)-vec2(2)*vec1(1);
			% If the cross product is negative, the angle is above pi and thus equal
			%	to 2pi-arcos(normalised scal). Otherwise, it's just arccos
			if (crossprod<0)
				coef=2*pi-acos(scal/(norm(vec1)*norm(vec2)));
			else
				coef=acos(scal/(norm(vec1)*norm(vec2)));
			end
	%coef

			sumcoefvert(j1)=sumcoefvert(j1)+coef;
			uvert(j1)=uvert(j1)+coef*uu(i);
		end
	end

	% Check if sumcoefvert is reasonable
	if (min(sumcoefvert)<1e-3 | max(sumcoefvert)>2*pi+1e-4)
		disp('wrong sumcoefvert');
	end
	%min(sumcoefvert)
	%max(sumcoefvert)
	% Then we apply the total ponderation
	uvert=uvert./sumcoefvert;

	% write file
	cd ensight/Un/;
	filename=strcat('Un',num2str(iter,'%05.f'));
	fid=fopen(filename,'w');
	fprintf(fid,'Un\npart\n1\ncoordinates\n');
	for i=1:ncell
		fprintf(fid,'%E\n',uu(i));
	end
	for j=1:nvert
		fprintf(fid,'%E\n',uvert(j));
	end;
	fclose(fid);
	
	cd ../..;

end





