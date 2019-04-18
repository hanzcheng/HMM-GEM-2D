%Computes the local flux matrix W_K such that
%	F_{K,sigma} = \sum_{sigma'} W^K_{sigma,sigma'}(u_K-u_{sigma'})
%
% xk is the x_K defined in the method, cg is the center of
% gravity of the cell - purely used to evaluate the diffusion inside the cell
function [W,G] = local_flux_matrix(vertex_loc,area,xk,cg);

% vertex_loc = (nb vertices + 1)x 2 array of vertex coordinates
%		vertex_loc(i,:) = coordinates of the i-th vertex of the cell K
%			The first vertex is repeated at the end of the matrix, so that
%			vertex_loc(i,:)-vertex_loc(i+1,:), for i=1,nb vertices, represents the
%			edges of the cell.
%		The vertices are oriented counter-clockwise

% Nb of edges in the cell
nbe = size(vertex_loc,1)-1;

G=sparse(2,nbe);
T=zeros(nbe,nbe);
X=sparse(nbe,2);
D=sparse(nbe,nbe);
W=sparse(nbe,nbe);

% midpoints of edges
xs=zeros(nbe,2);
xs([1:nbe],:)=(vertex_loc([1:nbe],:)+vertex_loc([2:nbe+1],:))/2;

% Outer normals: we start by rotating by -pi/2 the vectors vertex_loc(i+1,:)-vertex_loc(i,:)
% 	Because these are line vectors, we need to multiply on the right with the
%		transpose of the (-pi/2) rotation matrix
%		This rotation gives outer normal vectors of length the measure of the edges
N=zeros(nbe,2);
N([1:nbe],:)=(vertex_loc([2:nbe+1],:)-vertex_loc([1:nbe],:))*[0 -1;1 0];

% Computation of the length of the edges
msigma=sqrt(sum(N.^2,2));

% Normalisation of N
N=N./[msigma msigma];

% Matrix G_K
G=(-N.*[msigma msigma]./area)';

% Matrix X_K=xs-x_K
X=xs-ones(nbe,1)*xk;

% Diffusion in the cell
L=lambda(cg);

% Matrix D_K=diag(|sigma|/d_Ksigma Lambda n_Ks.n_Ks).
%		d_Ksigma=(xs-xk).n_Ksigma.
%% First option: simpler scaling by the trace of L
%D=diag(msigma./diag(X*N'))*trace(L);
%% Second option: the complete formula
% 		diag(N*(L*N')) is "(Lambda n_Ks.nKs)_{s edge of K}"
D=diag(msigma./diag(X*N').*diag(N*L*N'));

% Matrix T_K
T=-eye(nbe)-X*G; 

% Stab parameter
alpha=1;

% Matrix W
W=((area)*G'*L*G) +((alpha^2)*T'*D*T);



