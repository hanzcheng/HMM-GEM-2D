% Compute the L2 error on the function and the gradient
%
function [L2u,L2gradu]=compute_errors(u,area,center,vertex,cell_v,cell_e)

% u must be an ncell+nedge vector with cell unknowns 1:ncell, and edge
% unknowns ncell+1:ncell+edge.
ncell=size(area,1);

L2u=sqrt(sum(area.*(u(1:ncell)-ue(center)).^2))/sqrt(sum(area.*ue(center).^2));

L2gradu=0;
normgradu=0;
for i=1:ncell
	nbe=size(cell_v{i},2)-1;
	% Local vertices
	vertex_loc=vertex(cell_v{i},:);
	N=zeros(nbe,2);
	N([1:nbe],:)=(vertex_loc([2:nbe+1],:)-vertex_loc([1:nbe],:))*[0 -1;1 0];

	% Computation of the length of the edges
	msigma=sqrt(sum(N.^2,2));
	% Normalisation of N
	N=N./[msigma msigma];
	% Matrix G_K
	G=(-N.*[msigma msigma]./area(i))';
	% nabla_Ku
	nablaKu = -G*u(ncell+cell_e{i}([1:nbe]));

	[Du1,Du2]=derivu(center(i,:));
	L2gradu=L2gradu + area(i)*norm(nablaKu-[Du1;Du2])^2;
	normgradu=normgradu + area(i)*norm([Du1;Du2])^2;

end

L2gradu=sqrt(L2gradu)/sqrt(normgradu);



