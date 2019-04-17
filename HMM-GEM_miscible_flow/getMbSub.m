function sIntercept=getMbSub(subvertex,nbe)
mExt=(subvertex(2:nbe+1,2)-subvertex(1:nbe,2)) ./ (subvertex(2:nbe+1,1)-subvertex(1:nbe,1));
dExt= subvertex(1:nbe,2) - mExt .* subvertex(1:nbe,1);
dExt(abs(mExt)==Inf)=subvertex(abs(mExt)==Inf,1);
mExt(abs(mExt)==Inf)=Inf;
sIntercept=[mExt';dExt'];
end