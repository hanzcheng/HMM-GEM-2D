function intFlux = interiorFluxes_consistent(Fksigma,area,vertex,center)
s=size(Fksigma,1); %gives number of edges since Fksigma is edge fluxes
A=eye(s,s);
b=zeros(s,1);
aT=zeros(s,1);
    v1=vertex(1,:);
    v2=vertex(2,:);
    v3=center;
    aT(1)=0.5*det([v1-v3;v2-v3]);
for i=2:s
    A(i,i-1)=-1;
    v1=vertex(i,:);
    v2=vertex(i+1,:);
    v3=center;
    aT(i)=0.5*det([v1-v3;v2-v3]); %Area of Triangle i
    A(1,i-1)=(aT(i-1)+aT(i))/(2*area);
    b(i)=-Fksigma(i)+aT(i)/area*sum(Fksigma);
end 
A(1,s)=(aT(s)+aT(1))/(2*area);
intFlux=A\b;
end