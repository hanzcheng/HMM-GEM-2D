function KRDarcyU = kReconstruction(Fksigma,area,vertex,center,vType)
%returns the values of a,b1, and b2 where u=a[x1;x2]+[b1;b2]
%starts with the triangle corresponding to cell edge 1, 2, and so on
if vType=="KR"
intFlux = interiorFluxes(Fksigma,area,vertex,center); %corresponds to KR velocity
elseif vType=="C"
intFlux = interiorFluxes_consistent(Fksigma,area,vertex,center); %C-velocity
elseif vType=="A"
intFlux = interiorFluxes_aux(Fksigma,area,vertex,center); % A-velocity
end
s=size(intFlux,1); %gives number of edges since each intFlux corresponds to an edge
KRDarcyU=zeros(s,3);
for i=1:s
    A=zeros(s,s);
    b=zeros(s,1);
    v1=vertex(i,:);
    v2=vertex(i+1,:);
    v3=center;
    emid=0.5*[v1+v2; v2+v3; v3+v1]; %edge midpoints
    %to get the outward normal
    N = [v2-v1; v3-v2; v1-v3] * [0 -1;1 0];
    % Computation of the length of the edges
    msigma=sqrt(sum(N.^2,2));
    % Normalisation of N
    N=N./[msigma msigma];
    A=[sum(emid.*N,2) N];
    if i==1
        b=[Fksigma(i);intFlux(i);-intFlux(s)];
    else
        b=[Fksigma(i);intFlux(i);-intFlux(i-1)];
    end
    b=b./msigma;
    KRDarcyU(i,:)=(A\b)';
end
KRDarcyU(abs(KRDarcyU(:,:))<1e-10)=0;
KRDarcyU=KRDarcyU';
end