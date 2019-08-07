function [L2err,L1err,pL2err,pL1err,allLayers,nLayers] = compute_errors(solex,cCurr,volume,ncell)

L2err = (sqrt(volume'*(cCurr(1:ncell)-solex).^2));
pL2err = L2err/(sqrt(volume'*(solex).^2));
L1err = volume'*abs(solex-cCurr(1:ncell));
pL1err = L1err/(volume'*abs(solex));
allLayers = zeros(ncell,1);
nPct = 5;
for i=1:ncell
    if solex(i)~=0
        if abs(cCurr(i)-solex(i))>nPct/100*abs(solex(i))
            allLayers(i) = i;
        end
    else
        if abs(cCurr(i)-solex(i))>nPct/100
            allLayers(i) = i;
        end
    end
end
nz = find(allLayers>0);
nLayers = length(nz);