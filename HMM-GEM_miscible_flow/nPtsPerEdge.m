function [pPerEdge,mReg] = nPtsPerEdge(tstep,diam,area,ncell)
pPerEdge=zeros(ncell,1);
mReg1=(diam .* diam) ./ area;
for i=1:ncell
    pPerEdge(i)=log(mReg1(i))/log(2);
    if abs(1-abs(pPerEdge(i)-ceil(pPerEdge(i))))<1e-10
        pPerEdge(i) = ceil(pPerEdge(i))-1;
    else
        pPerEdge(i) = ceil(pPerEdge(i));
    end
    %% can comment out the latter lines for mKreg points
        pPerEdge(i) = 2^(ceil(2*tstep/diam(i)))*pPerEdge(i)+1;
        if mod(pPerEdge(i),2)==0
            pPerEdge(i) = pPerEdge(i)+1; % odd no. of points to make sure midpt is included in tracking
        end
end
mReg = max(mReg1);
end