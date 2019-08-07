function solex = compute_refSol(v_exact,tFin,ncell,cell_v,vertex,center)

%% note: by choosing snQuadPts = 2^n, it is as if we have the sol. from refining the mesh by n times
snQuadPts = 4; %sqrt of the no. of quadrature points 
[cB,mainCell] = compute_qPoints(ncell,cell_v,vertex,snQuadPts);
pointsToTrack = cB;
%edge midpoints

tstep=1e-1;
nbSteps=ceil(tFin/tstep);

solini=initCond(center);
solex = zeros(size(solini));

solex_ref = initCond(track_exact(pointsToTrack,nbSteps*tstep,v_exact)); % "exact" sol. on a refined grid
solex_ref(solex_ref==1)=1;
for i=1:ncell
    cellsNow = find(mainCell==i);
    nCellsNow = length(cellsNow);
    solex(i) = sum(solex_ref(cellsNow))/nCellsNow;
end

end