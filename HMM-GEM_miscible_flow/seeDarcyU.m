function seeDarcyU(sub_cell_e,sub_cell_n,sub_cell_v,sub_vertex,sub_N,sub_ncell,sub_nvert,sub_triInvVert,sIntercept,cornerPoints,KRDarcyU)
% cellsStart=[];
% for i=1:size(icell,2)
% cellsStart=[cellsStart cell_n{icell(i)}(cell_n{icell(i)}~=0) icell(i)];
% end
% nStarting=size(cellsStart,2);
% sLoc=[];
% for i=1:nStarting
%     sLoc=[sLoc find(mainCell==cellsStart(i))];
% end
myTri=1:sub_ncell;
vertices=cell(1,sub_ncell);
darcyUHere=KRDarcyU(:,myTri);
cv=zeros(sub_ncell,2);
u=zeros(sub_ncell,1);
v=zeros(sub_ncell,1);
x=zeros(sub_ncell,1);
y=zeros(sub_ncell,1);
for i=1:sub_ncell
    vertices{i}=sub_vertex(sub_cell_v{myTri(i)},:);
    cv(i,:)=sum(vertices{i}(1:3,:))/3;
    a=darcyUHere(1,i);
    b1=darcyUHere(2,i);
    b2=darcyUHere(3,i);
    u(i)=a*cv(i,1)+b1;
    v(i)=a*cv(i,2)+b2;
    x(i)=cv(i,1);
    y(i)=cv(i,2);
end
quiver(x,y,u,v)
locStart1=intersect(find(x<1000),find(x>900));
locStart2=intersect(find(y<1000),find(y>900));
locStart=intersect(locStart1,locStart2);
startx=x(locStart);
starty=y(locStart);
pointsToTrack=[startx starty];
triInvI=locStart;
tstep=36;
nbStream=100;
% nbStream=70;
nbPts=size(locStart,1);
xy=cell(nbPts,1);
for i=1:nbPts
    xy{i}=[xy{i} ; pointsToTrack(i,:)];
end
for j=1:nbStream
    [~,newPointsF,~,nextTriF]=complete_tracking_efficient(sub_cell_e,sub_cell_n,sub_cell_v,sub_vertex,sub_N,sub_nvert,sub_triInvVert,sIntercept,tstep,cornerPoints,pointsToTrack,1:nbPts,1:nbPts,triInvI,KRDarcyU);
    pointsToTrack=newPointsF;
    triInvI=nextTriF;
    for i=1:nbPts
        if pointsToTrack(i,1)>0 && pointsToTrack(i,2)>0 && pointsToTrack(i,1)<1000 && pointsToTrack(i,2)<1000
            xy{i}=[xy{i} ; pointsToTrack(i,:)];
        end
    end
end
streamline(xy)

%     title('Reconstructed Darcy Velocity at initial time t=0 (Kershaw mesh)')
xlabel('x')
ylabel('y')

end