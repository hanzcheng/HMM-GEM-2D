function [newPointsB,newPointsF,nextTriB,nextTriF]=complete_tracking_efficient(sub_cell_e,sub_cell_n,sub_cell_v,sub_vertex,sub_N,sub_nvert,sub_triInvVert,sIntercept,tstep,cornerPoints,pointsToTrack,pNearInj,pNearProd,triInvI,KRDarcyU)
bTracking=0;
fTracking=1;

%tracking is 1 if we will do forward tracking, and 0 for backward tracking
time_leftB=tstep*ones(length(pNearProd),1); %time_end is t^{n+1}
time_leftF=tstep*ones(length(pNearInj),1);
locCorner=checkPoints(cornerPoints,pointsToTrack);
%% initialize triangles for tracking
currTriB=zeros(length(pNearProd),1);
% tic;
%% backward
% determine the triangles for the initial set of points to track
currTriB=detTri(triInvI(pNearProd,:),currTriB,sub_cell_n,sub_cell_v,sub_vertex,sub_N,KRDarcyU,pointsToTrack(pNearProd,:),bTracking,locCorner);
[time_To_ExitB,nextTriB]=tExit(currTriB,sub_cell_e,sub_cell_n,sIntercept,KRDarcyU,time_leftB,pointsToTrack(pNearProd,:),bTracking);

%% backward tracking
newPointsB=pointsToTrack(pNearProd,:);

locToBTrack=find(time_leftB>0);%location(s) of the point(s) that are still to be tracked since tstep has not yet been fully utilized

while ~all(time_leftB<=0)
    triInvNowB=zeros(length(pointsToTrack),sub_nvert);
    %track
    newPointsB(locToBTrack,:)=charTracking1(newPointsB(locToBTrack,:),time_To_ExitB(locToBTrack),currTriB(locToBTrack),KRDarcyU,bTracking);
    %determine next triangle
    triInvNowB(locToBTrack,1)=nextTriB(locToBTrack);
    %now we check if the newPoints coincide with one of the vertices
    posSame=checkPoints(newPointsB(locToBTrack,:),sub_vertex);
    newPointsB(locToBTrack(posSame~=0),:)=sub_vertex(posSame(posSame~=0),:);
    locCorner=checkPoints(cornerPoints,newPointsB(locToBTrack,:));
    triInvNowB(locToBTrack(posSame~=0),:)=sub_triInvVert(posSame(posSame~=0),:);
    %compute time left
    time_leftB(locToBTrack)=time_leftB(locToBTrack)-time_To_ExitB(locToBTrack);
    %check which point(s) has yet to be tracked
    locToBTrack=find(time_leftB>0);
    %update the triangle
    currTriB(locToBTrack)=detTri(triInvNowB(locToBTrack,:),currTriB(locToBTrack),sub_cell_n,sub_cell_v,sub_vertex,sub_N,KRDarcyU,newPointsB(locToBTrack,:),bTracking,locCorner);
    nextTriB(time_leftB<=0)=currTriB(time_leftB<=0);
    %compute time of exit
    [time_To_ExitB(locToBTrack),nextTriB(locToBTrack)]=tExit(currTriB(locToBTrack,:),sub_cell_e,sub_cell_n,sIntercept,KRDarcyU,time_leftB(locToBTrack),newPointsB(locToBTrack,:),bTracking);
    
end
posSame=checkPoints(newPointsB(locToBTrack,:),sub_vertex);
newPointsB(locToBTrack(posSame~=0),:)=sub_vertex(posSame(posSame~=0),:);


%% forward tracking
currTriF=zeros(length(pNearInj),1);

%% forward
% determine the triangles for the initial set of points to track
currTriF=detTri(triInvI(pNearInj,:),currTriF,sub_cell_n,sub_cell_v,sub_vertex,sub_N,KRDarcyU,pointsToTrack(pNearInj,:),fTracking,locCorner);
[time_To_ExitF,nextTriF]=tExit(currTriF,sub_cell_e,sub_cell_n,sIntercept,KRDarcyU,time_leftF,pointsToTrack(pNearInj,:),fTracking);


newPointsF=pointsToTrack(pNearInj,:);

locToFTrack=find(time_leftF>0);%location(s) of the point(s) that are still to be tracked since tstep has not yet been fully utilized

while ~all(time_leftF<=0)
    triInvNowF=zeros(length(pointsToTrack),sub_nvert);
    newPointsF(locToFTrack,:)=charTracking1(newPointsF(locToFTrack,:),time_To_ExitF(locToFTrack),currTriF(locToFTrack),KRDarcyU,fTracking);
    triInvNowF(locToFTrack,1)=nextTriF(locToFTrack);
    %check if the newPoints coincide with one of the vertices
    posSameF=checkPoints(newPointsF(locToFTrack,:),sub_vertex);
    newPointsF(locToFTrack(posSameF~=0),:)=sub_vertex(posSameF(posSameF~=0),:);
    locCorner=checkPoints(cornerPoints,newPointsF(locToFTrack,:));
    triInvNowF(locToFTrack(posSameF~=0),:)=sub_triInvVert(posSameF(posSameF~=0),:);
    time_leftF(locToFTrack)=time_leftF(locToFTrack)-time_To_ExitF(locToFTrack);
    locToFTrack=find(time_leftF>0);
    currTriF(locToFTrack)=detTri(triInvNowF(locToFTrack,:),currTriF(locToFTrack),sub_cell_n,sub_cell_v,sub_vertex,sub_N,KRDarcyU,newPointsF(locToFTrack,:),fTracking,locCorner);
    
    nextTriF(time_leftF<=0)=currTriF(time_leftF<=0);
    [time_To_ExitF(locToFTrack),nextTriF(locToFTrack)]=tExit(currTriF(locToFTrack,:),sub_cell_e,sub_cell_n,sIntercept,KRDarcyU,time_leftF(locToFTrack),newPointsF(locToFTrack,:),fTracking);
    
end
posSameF=checkPoints(newPointsF(locToFTrack,:),sub_vertex);
newPointsF(locToFTrack(posSameF~=0),:)=sub_vertex(posSameF(posSameF~=0),:);

% tTrack=toc;