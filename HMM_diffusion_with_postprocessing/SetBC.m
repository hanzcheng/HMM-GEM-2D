% Set the boundary conditions on the cell i
%
%		bccase=0: pure Dirichlet
%		bccase=1: pure Neumann (NOT IMPLEMENTED YET)
%		bccase>1: mixed
%
%	The boundary condition on the edge cell_e{i}(j) is imposed by
%		setting cell_n{i}(j) at 0 for Dirichlet BC, -1 for Neumann BC
%	We actually modify cell_n{i} where needed
%
function bccell=SetBC(cell_n,vertex)

global bccase;

% Nb of edges in the cell
nbe = size(vertex,1)-1;

bccell=cell_n;

if (bccase==0)
	% Pure Dirichlet, no change in cell_n
elseif (bccase==1)
    % Pure Neumann
    for j=1:nbe
        if(bccell(j)==0)
            bccell(j)=-1;
        end
    end
elseif (bccase==2)
	% Neumann at x=1, Dirichlet elsewhere
	for j=1:nbe
		% We only change Neumann edges
		if ( (vertex(j,1)+vertex(j+1,1))/2 > 1-1e-5)
			bccell(j)=-1;
		end
	end
end



