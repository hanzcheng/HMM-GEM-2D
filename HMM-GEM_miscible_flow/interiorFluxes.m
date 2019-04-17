function intFlux = interiorFluxes(Fksigma,area,vertex,center)
s=size(Fksigma,1); %gives number of edges since Fksigma is edge fluxes
A=eye(s,s);
b=zeros(s,1);
for i=2:s
    A(i,i-1)=-1;
    v1=vertex(i,:);
    v2=vertex(i+1,:);
    v3=center;
     aT=0.5*det([v1-v3;v2-v3]); %Area of Triangle i
    b(i)=-Fksigma(i)+aT/area*sum(Fksigma);
%     if(abs(b(i))<1e-10)
%         b(i)=0;
%     end
end
A=A(2:s,:);
b=b(2:s);
intFlux=A'*((A*A')\b);
end