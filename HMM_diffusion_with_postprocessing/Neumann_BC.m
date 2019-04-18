function gv=Neumann_BC(v1,v2,L)
cv = 0.5 * (v1+v2); %gets the center of each edge
%to get the outward normal
N = (v2-v1) * [0 -1;1 0];

% Computation of the length of the edges
msigma=sqrt(sum(N.^2,2));
% Normalisation of N
N=N./[msigma msigma];
ne=size(cv,1);
%gives the number of edges
gv =zeros(ne,1);
x=N(:,1);
y=N(:,2);
for i=1:ne
    
	[gx,gy]=derivu(cv(i,:));
	Lg=L*[gx;gy];
	Lgx=Lg(1);
	Lgy=Lg(2);
	gv(i) =(Lgx*x(i)+Lgy*y(i))*msigma(i);

end
%gv=-(2*y*(y-1)+2*x*(x-1)-(x-x.^2).*(y-y.^2));
%gv=2*pi^2*sin(pi*x)*sin(pi*y);
%gv=pi*x*sin(pi*x)*sin((pi*y)/2) - 2*pi*y*cos(pi*x)*cos((pi*y)/2) + (5*pi^2*x*y*cos((pi*y)/2)*sin(pi*x))/4;
