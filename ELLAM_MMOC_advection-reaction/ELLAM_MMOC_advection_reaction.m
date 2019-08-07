%% Characteristic-based schemes for the advection-reaction problem
% for solenoidal fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_t + alpha*c + div (cV) = f
% alpha is the reaction coefficient
% V is the velocity with div(V)=0
% with source term f and initial condition c(x,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% -------------------------------------------------------------------
% Setup and load mesh
mesh = {'mesh2_3.mat'};       % The mesh file in MATLAB data format

loadmesh=strcat('load ../../matlab_meshes/',mesh{1});
eval(loadmesh);
cg = gravity_centers(ncell,cell_v,vertex,area);
% -------------------------------------------------------------------
% Setup problem parameters
% reaction coefficient
alfa = 0;
% Source term
source = @(x,y) 0;
global testcase
testcase = 4;
% Exact velocity
if testcase <= 2
    v_exact = @(x,y) [ 0.0625 *ones(size(x))  0 * ones(size(y))];
elseif testcase >= 3 
    v_exact = @(x,y) [(1-2*y) .* x .* (1-x), -(1-2*x) .* y .* (1-y)];
end
% Initial condition

tFin = 8; %final time
tstep = 0.8; %time step
solex = compute_refSol(v_exact,tFin,ncell,cell_v,vertex,center);

cPrev = zeros(ncell+nedge,1);
for i=1:ncell
    cPrev(i) = initCond(cg(i,:));
end

%% choice of numerical scheme and time step
numScheme = "ELLAM"; %actually since we look at solenoidal fields
% numScheme = "MMOC"; %ELLAM and MMOC are equivalent


nbsteps = ceil(tFin/tstep);
times=[0:tstep:tFin];

% determine no. of points to track per edge
pPerEdge = nPtsPerEdge(tstep,diam,area,ncell);
MpPerEdge = max(pPerEdge);

pOnEdge=zeros(nedge,1); %pOnEdge(i) gives the number of points to be tracked along edge i
for cell_i=1:ncell
    nbe=length(cell_e{cell_i});
    for j=1:nbe
        pOnEdge(cell_e{cell_i}(j))=max(pOnEdge(cell_e{cell_i}(j)),pPerEdge(cell_i));
    end
end

[xEdges,yEdges]= edge_points(ncell,cell_e,cell_v,vertex,pPerEdge);
%for documentation on xEdges, see edge_points
nxEdges=cell(1,ncell);
nyEdges=cell(1,ncell);
pointsToTrack=zeros(ncell*nedge*(MpPerEdge+2)+nvert,2);
ctr=1;

%here, we modify the number of points on each edge so that edges shared between 2 cells have matching pts along edges
for cell_i=1:ncell
    currNeigh=cell_n{cell_i};
    nbe=length(cell_e{cell_i});
    nxEdges{cell_i}=cell(1,nbe);
    for j=1:nbe
        if size(xEdges{cell_i}(:,j),1)<pOnEdge(cell_e{cell_i}(j))+2
            copyNeigh=cell_n{cell_i}(j); %neighbor to copy from
            copyEdge=find(cell_n{copyNeigh}==cell_i); %edge to copy from
            nxEdges{cell_i}{j}=flip(xEdges{copyNeigh}(:,copyEdge));
            nyEdges{cell_i}{j}=flip(yEdges{copyNeigh}(:,copyEdge));
        else
            nxEdges{cell_i}{j}=xEdges{cell_i}(:,j);
            nyEdges{cell_i}{j}=yEdges{cell_i}(:,j);
        end
    end
    
end
for cell_i=1:ncell
    nbe=length(cell_e{cell_i});
    for j=1:nbe
        pointsToTrack(ctr:ctr+pOnEdge(cell_e{cell_i}(j))-1+2,:)=[nxEdges{cell_i}{j} nyEdges{cell_i}{j}];
        ctr=ctr+pOnEdge(cell_e{cell_i}(j))+2;
    end
end

pointsToTrack=uniquetol(pointsToTrack,'ByRows',true);
xm=min(pointsToTrack(:,1));
xM=max(pointsToTrack(:,1));
ym=min(pointsToTrack(:,2));
yM=max(pointsToTrack(:,2));
cornerPoints=[xm ym; xm yM; xM ym; xM yM];

cellInvI=zeros(size(pointsToTrack,1),ncell); %initial cells involved per point
%cell_ep(i,:) gives the location/number of the points(to be tracked)
%that are located along cell i of the mesh,i.e. which points are in
%cell i?
cell_ep=whichPoint_Cell(ncell,cell_e,nxEdges,nyEdges,pointsToTrack,MpPerEdge+2);
for cell_i=1:ncell
    for s=1:length(cell_ep(cell_i,:))
        if cell_ep(cell_i,s)~=0
            posZ=find(cellInvI(cell_ep(cell_i,s),:)==0);
            cellInvI(cell_ep(cell_i,s),posZ(1))=cell_i;
        end
    end
end


% characteristic tracking
if numScheme == "ELLAM"
    [newPoints,tTrack]=complete_tracking_exact_vel(tstep,pointsToTrack,v_exact,0);
elseif numScheme == "MMOC"
    [newPoints,tTrack]=complete_tracking_exact_vel(tstep,pointsToTrack,v_exact,1);
end
% area of intersection bet tracked and residing cells
[AreaIn,AreaT,tInt]=compArea_cell(ncell,pointsToTrack, newPoints,cell_ep);
% local volume adj for local mass conservation
tic;
if testcase>=3
[AreaIn,pctErr,aErr,nAdj,aErrTbInit] = adjust_volumes_for_local_mass_divFree(vertex,cell_v,AreaIn,area,ncell,cell_e,cell_n,newPoints,cellInvI,cell_ep,pOnEdge,v_exact);
end
tAdj=toc;
tTotal = tAdj+tInt+tTrack;
tTrack = tTrack/tTotal;
tAdj = tAdj/tTotal;
tInt = tInt/tTotal;

write_solution_ensight(0,0,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
%process the initial condition for plotting by 
% setting the value at the edges to be the average of the value of 2
% adjacent cells 
    for i=1:ncell
        nbe = length(cell_e{i});
        for j=1:nbe
            if cPrev(ncell+cell_e{i}(j))==0
                cPrev(ncell+cell_e{i}(j)) = cPrev(i);
            else
                cPrev(ncell+cell_e{i}(j)) =0.5*( cPrev(ncell+cell_e{i}(j))+cPrev(i));
            end
        end
    end
write_solution_ensight(1,cPrev,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);

%measure of mass conservation
massCons = zeros(1,nbsteps);
%% solve the advection PDE
for m=1:nbsteps
    A = zeros(ncell,1);
    b = zeros(ncell,1);
    
    % assemble LHS and source term
    for cell_i=1:ncell
        A(cell_i) = area(cell_i);
        if numScheme == "ELLAM"
            intC = AreaIn(:,cell_i)'*cPrev(1:ncell);
        elseif numScheme == "MMOC"
            intC = AreaIn(cell_i,:)*cPrev(1:ncell);
        end
        b(cell_i) = intC +tstep* (source(cg(cell_i,1),cg(cell_i,2)) - alfa*cPrev(cell_i));
    end

    % solve for c
    cCurr = b ./ A;
    
    %measure of global mass conservation bet. current time and prev. time
    massCons(m) = sum(area.*cPrev(1:ncell))-sum(area.*cCurr(1:ncell));
    
    cPrev = cCurr;
    c = zeros(ncell+nedge,1);
    c(1:ncell) = cPrev;
    % set the value at the edges to be the average of the value of 2
    % adjacent cells (to be used for plotting)
    for i=1:ncell
        nbe = length(cell_e{i});
        for j=1:nbe
            if c(ncell+cell_e{i}(j))==0
                c(ncell+cell_e{i}(j)) = cCurr(i);
            else
                c(ncell+cell_e{i}(j)) =0.5*( c(ncell+cell_e{i}(j))+cCurr(i));
            end
        end
    end
    c(abs(c)<1e-12)=0;
    
    write_solution_ensight(1,c,m,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
end
[L2err,L1err,pL2err,pL1err,allLayers,nLayers] = compute_errors(solex,cCurr,area,ncell);
    L2err
    L1err
    pL2err
    pL1err

if numScheme == "ELLAM"
    write_solution_vtk(c,'solution_ELLAM_rect',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
elseif numScheme == "MMOC"
    write_solution_vtk(c,'solution_MMOC_rect',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
end
    