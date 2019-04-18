function [newPoints,tTrack]=complete_tracking_exact_vel(tstep,pointsToTrack,dv_exact,tTracking)
%% function that returns the location of the points after completely tracking tstep
% inputs: 
%       : tstep -- time step to be taken
%       : pointsToTrack -- points to be tracked
%       : dv_exact -- exact Darcy velocity
%       : tTracking -- 0 is backward tracking, 1 is forward tracking
% outputs: newPoints -- location of the points that are tracked
%        : tTrack -- time it takes to perform the complete tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subTsteps =100; % since using first order Euler, accuracy/precision of the 
               % tracking is made by taking sub-timesteps, in this case 100 
time_leftB=tstep*ones(length(pointsToTrack),1); %time_end is t^{n+1}
time_To_ExitB = time_leftB/subTsteps;
tic;
%% backward
newPoints=pointsToTrack;
locToBTrack=find(time_leftB>0); %location(s) of the point(s) that are 
           %still to be tracked since tstep has not yet been fully utilized
           
while ~all(time_leftB<=0)
        ctrTstep = 0; %counter for no. of time steps taken
        while ctrTstep<subTsteps
            newPoints(locToBTrack,:)=FO_Euler_exactVel(newPoints(locToBTrack,:),time_To_ExitB(locToBTrack),dv_exact,tTracking);
            time_leftB(locToBTrack) = time_leftB(locToBTrack)-time_To_ExitB(locToBTrack);
            locToBTrack=find(time_leftB>0);
            ctrTstep = ctrTstep+1;
        end
 
end
tTrack=toc;