function fv=source_term(z)
[ux,uy,H]=derivu(z);
[L,dL]=lambda(z);
Lh=L .* H;
fv=-(sum(sum(Lh))+dL(1,1)*ux+dL(1,2)*uy+dL(2,1)*ux+dL(2,2)*uy);
