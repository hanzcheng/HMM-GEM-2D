% Track according to minus velocity

function xT=track_exact(x0,T,vel);

dt=1e-3;
Ndt=T/dt;


xT=x0;
for idt=1:Ndt
    xT = xT - dt*vel(xT(:,1),xT(:,2));
%     for j=1:length(xT)
% 	xT(j,:) = xT(j,:) - dt * vel(xT(j,1),xT(j,2))';
%     end
end;



