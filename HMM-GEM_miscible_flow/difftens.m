function D=difftens(u,Phi,Phidm,Phidl,Phidt)
%gives the diffusion tensor D
u1=u(1);
u2=u(2);
uNorm=sqrt(u1^2+u2^2);
proj=1/uNorm^2*[u1*u1, u1*u2; u2*u1, u2*u2];
D=(Phidm*eye(2,2)+uNorm*(Phidl*proj+Phidt*(eye(2,2)-proj)));
end