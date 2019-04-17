%   Hanz Martin Cheng 
%   hanz.cheng@monash.edu
%   Monash University, Australia
%   
%   HMM-MMOC scheme for the miscible flow model on generic meshes
%
%   Unknowns: pressure and concentration.
%   Auxliary unknown: darcy velocity
%
%   Test cases from the paper by Wang et al,
%   SIAM, J. Sci. Comput.
%
%   Acknowledgement: This research was supported by the 
%   Australian Government through the Australian Research Council's 
%   Discovery Projects funding scheme (project number DP170100605)
%
%   This Source Code Form is subject to the terms of the GNU Lesser General
%   Public License v3.0. If a copy of the LGPL was not distributed with 
%   this file, you can obtain one at https://www.gnu.org/licenses/lgpl-3.0
%  
%   If you use this code or parts of it for scientific publications, you
%   are required to cite it as following:
%   
%   GEM component and local volume adjustments:
%   A combined GDM--ELLAM--MMOC scheme for advection dominated PDEs
%   H.M. Cheng, J. Droniou and K.-N. Le
%   ArXiv e-prints
%   http://adsabs.harvard.edu/abs/2018arXiv180505585C
%
%   ELLAM component:
%   An HMM--ELLAM scheme on generic polygonal meshes for miscible
%   incompressible flows in porous media. 
%   H.M. Cheng and J. Droniou
%   Journal of Petroleum Science and Engineering
%   DOI: 10.1016/j.petrol.2018.08.062
%
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMAIN AND MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
mesh={'mesh2_3.mat'}; %rectangular mesh
%    mesh={'mesh4_1.mat'}; %Kershaw mesh
% mesh={'hexa1_2.mat'}; %hexahedral mesh
% mesh={'non_conforming_3.mat'};
loadmesh=strcat('load ../../../matlab_meshes/',mesh{1});
eval(loadmesh);
cg = gravity_centers(ncell,cell_v,vertex,area);

diam=diam*1000;
area=area*1000^2;
vertex=vertex*1000;
center=center*1000;
cg=cg*1000;

%% Diffusion choice
diffType = 1; %standard diffusion tensor
% diffType = 2; %vanishing diffusion (to help reduce artificial fingers at bdry)

%% time step
tstep=36;

%% create submesh consisting of triangles (to be used later for RT0)
[sub_ncell,sub_nvert,sub_nedge,sub_vertex,sub_cell_v,sub_cell_n,sub_cell_e,sub_N,sub_triInvVert,mainCell,subCellsIn]=create_submesh(ncell,nvert,vertex,cg(:,:),cell_v,cell_e);
vType = "A"; %type of velocity reconstruction: can choose A, C, or KR
% vType = "C";
% vType = "KR";

sIntercept=zeros(2,sub_nedge);
% to obtain the slope and intercept of the segment that determines the
% edges of subcell i, we simply apply sIntercept(:,sub_cell_e{i}) 1st row is
% the slope, while 2nd row is the intercept
for i=1:sub_ncell
    sIntercept(:,sub_cell_e{i})=getMbSub(sub_vertex(sub_cell_v{i},:),3);
end

%% determine number of points to track along the edge of each polygon
[pPerEdge,mReg] = nPtsPerEdge(tstep,diam,area,ncell);
% mReg gives a measure of the regularity of the mesh
%% alternatively, if we want to track a fixed no. of pts per edge, can use
%N = 3;
% for i=1:ncell
%     pPerEdge(i)=N;
% end
MpPerEdge=max(pPerEdge);

qinj=30;  %injection at corner (1000,1000)
qprod=30; %production at corner(0,0)
cinj=1.00;

%functions associated to sources

qm=zeros(ncell+nedge,1);
qp=zeros(ncell+nedge,1);
chat=zeros(ncell+nedge,1); %chat = c "hat"
xprod=0;
yprod=0;    % x and y coordinates of the production point
xinj=1000;
yinj=1000;  % x and y coordinates of the injection point
thres=min(diam)/(10*(MpPerEdge+2)); %threshold is 10% of the mesh size

%% locate the injection and production cells
[icell,pcell,iCellCount,pCellCount]= locateWells(xprod,yprod,xinj,yinj,vertex,cell_v,ncell,thres);

qp(icell)=qinj * area(icell) / sum(area(icell));
qm(pcell)=qprod * area(icell) / sum(area(icell));
q=qp-qm;
chat(icell)=cinj;
pPerEdge(icell)=4*MpPerEdge+1; %need to track more pts at icell and pcell
pPerEdge(pcell)=4*MpPerEdge+1;
%%
pOnEdge=zeros(nedge,1); %pOnEdge(i) gives the number of points to be tracked along edge i
for i=1:ncell
    nbe=length(cell_e{i});
    for j=1:nbe
        pOnEdge(cell_e{i}(j))=max(pOnEdge(cell_e{i}(j)),pPerEdge(i));
    end
end

[xEdges,yEdges]= edge_points(ncell,cell_e,cell_v,vertex,pPerEdge);
%for documentation on xEdges, see edge_points
nxEdges=cell(1,ncell);
nyEdges=cell(1,ncell);
pointsToTrack=zeros(ncell*nedge*(MpPerEdge+2)+nvert,2);
ctr=1;
mnbe=0;

%here, we modify the number of points on each edge so that edges shared between 2 cells have matching pts along edges
for i=1:ncell
    currNeigh=cell_n{i};
    nbe=length(cell_e{i});
    nxEdges{i}=cell(1,nbe);
    for j=1:nbe
        if size(xEdges{i}(:,j),1)<pOnEdge(cell_e{i}(j))+2
            copyNeigh=cell_n{i}(j); %neighbor to copy from
            copyEdge=find(cell_n{copyNeigh}==i); %edge to copy from
            nxEdges{i}{j}=flip(xEdges{copyNeigh}(:,copyEdge));
            nyEdges{i}{j}=flip(yEdges{copyNeigh}(:,copyEdge));
        else
            nxEdges{i}{j}=xEdges{i}(:,j);
            nyEdges{i}{j}=yEdges{i}(:,j);
        end
    end
    
end
for i=1:ncell
    nbe=length(cell_e{i});
    mnbe=max(mnbe,nbe); %to get the maximum number of edges in the submesh
    for j=1:nbe
        pointsToTrack(ctr:ctr+pOnEdge(cell_e{i}(j))-1+2,:)=[nxEdges{i}{j} nyEdges{i}{j}];
        ctr=ctr+pOnEdge(cell_e{i}(j))+2;
    end
end

pointsToTrack=uniquetol(pointsToTrack,'ByRows',true);
xm=min(pointsToTrack(:,1));
xM=max(pointsToTrack(:,1));
ym=min(pointsToTrack(:,2));
yM=max(pointsToTrack(:,2));
cornerPoints=[xm ym; xm yM; xM ym; xM yM];

cell_ep=cell(1,ncell);
triInvI=zeros(size(pointsToTrack,1),sub_ncell); %initial triangles involved per point
%sub_cell_ep(n,:) gives the location/number of the points(to be tracked)
%that are located along triangle n of the submesh
sub_cell_ep=whichPoint(ncell,sub_ncell,cell_e,nxEdges,nyEdges,pointsToTrack,MpPerEdge+2);
for i=1:sub_ncell
    for s=1:length(sub_cell_ep(i,:))
        if sub_cell_ep(i,s)~=0
            posZ=find(triInvI(sub_cell_ep(i,s),:)==0);
            triInvI(sub_cell_ep(i,s),posZ(1))=i;
            
            
            posZ=find(triInvI(sub_cell_ep(i,s),:)==0);
            mCells=mainCell(sub_cell_n{i}(sub_cell_n{i}~=0));
            uTri=sub_cell_n{i}(mCells~=mainCell(i)); %undetected sub triangle
            if ~isempty(uTri)
                triInvI(sub_cell_ep(i,s),posZ(1))=uTri;
            end
            
        end
    end
end

%to form the approximation of polygon i
%obtained after doing the tracking, we do the following
%1. get length of subCellsIn{i}
%2. loop over subCellsIn{i}(j) and take
%newPoints(sub_cell_ep(subCellsIn{i}(j),:))
%detect boundary edges for main mesh
bdI=cell(1,ncell);
for i=1:ncell
    bdI{i}=find(cell_n{i}==0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

castest=2;

switch castest
    case 1
        %viscosity with mu0=1 and M=41
        mu = @(c) ones(size(c));
        %porosity
        Phi=0.1;
        
        %permeability
        kappa=@(x,y) 80.*(y>=0);
        
        %molecular diffusion
        Phidm=1.0;
        
        %dispersions (longitudinal and transverse)
        Phidl=0;
        Phidt=0;
        
    case 2
        mu = @(c) (1-c+41^(1/4)*c).^(-4);
        
        Phi=0.1;
        
        kappa=@(x,y) 80.*(y>=0);
        
        Phidm=0.0;
        
        Phidl=5.0;
        Phidt=0.5;
        
    case 3
        mu = @(c) (1-c+41^(1/4)*c).^(-4);
        
        Phi=0.1;
        
        kappa = @(x,y) 800*(y<500)+2*(y>=500);
        
        Phidm=0.0;
        
        Phidl=5.0;
        Phidt=0.5;
        
    case 4
        mu = @(c) (1-c+41^(1/4)*c).^(-4);
        
        Phi=0.1;
        
        kappa = @(x,y) 80-55*(x>150).*(x<550).*(y>150).*(y<550);
        
        Phidm=0.0;
        
        Phidl=5.0;
        Phidt=0.5;
        
        
    case 5
        mu = @(c) (1-c+41^(1/4)*c).^(-4);
        
        Phi=0.1;
        
        kappa = @(x,y) 80-60*(sin(5*pi*x/1000)<0).*(sin(5*pi*y/1000)<0);
        
        Phidm=0.0;
        
        Phidl=5.0;
        Phidt=0.5;
        
    case 6
        mu = @(c) (1-c+41^(1/4)*c).^(-4);
        
        Phi=0.1;
        
        kappa = @(x,y) 80-60*(x>200).*(x<400).*(y>200).*(y<400)-60*(x>200).*(x<400).*(y>600).*(y<800)-60*(x>600).*(x<800).*(y>600).*(y<800)-60*(x>600).*(x<800).*(y>200).*(y<400);
        
        Phidm=0.0;
        
        Phidl=5.0;
        Phidt=0.5;
end

%% Pre-determine MMOC regions by approximating which cells are nearer to injection
allCells=1:ncell;
nearInjCell=allCells;
currDist1=1000*sqrt(2)*ones(1,ncell);
currDist2=1000*sqrt(2)*ones(1,ncell);
nearInjCellIni=zeros(1,ncell);
for i=1:ncell
    for j=1:iCellCount
        currDist1(i)=min(currDist1(i),sqrt(sum(center(i,:)-center(icell(j),:)).^2));
    end
end
for i=1:ncell
    for j=1:pCellCount
        currDist2(i)=min(currDist2(i),sqrt(sum(center(i,:)-center(pcell(j),:)).^2));
    end
    
    if currDist1(i)<currDist2(i) %nearer to the injection well, trace forward
        nearInjCellIni(i)=i;
    end
end
nearProdCell=allCells-nearInjCell;
nearInjCellIni(nearInjCellIni==0)=[];
nearProdCell(nearProdCell==0)=[];
%%
TriNearInj=zeros(1,mnbe*length(nearInjCellIni));%triangles near the injection cells
ctr=1;
for i=1:length(nearInjCellIni)
    nbe=length(cell_e{nearInjCellIni(i)});
    TriNearInj(ctr:ctr+nbe-1)=subCellsIn{nearInjCellIni(i)}; %to get the triangles involved with the injection well(s)
    ctr=ctr+nbe;
end
TriNearInj(TriNearInj==0)=[];
pNearInjIni=sub_cell_ep(TriNearInj,:);
pNearInjIni(pNearInjIni==0)=[];
pNearInjIni=unique(pNearInjIni(:));
pNearInj=1:length(pointsToTrack);
%%
TriNearProd=zeros(1,mnbe*length(nearProdCell));%triangles near the production cells
ctr=1;
for i=1:length(nearProdCell)
    nbe=length(cell_e{nearProdCell(i)});
    TriNearProd(ctr:ctr+nbe-1)=subCellsIn{nearProdCell(i)}; %to get the triangles involved with the production well(s)
    ctr=ctr+nbe;
end
TriNearProd(TriNearProd==0)=[];
pNearProd=sub_cell_ep(TriNearProd,:);
pNearProd(pNearProd==0)=[];
pNearProd=unique(pNearProd(:));

pNearProdIni=1:length(pointsToTrack);
%%
tstep=36; %time step
nbyrs=10; %number of years
nbiter=ceil(nbyrs*360/tstep); %number of iterations
alffa=sum(qp(icell))*tstep/(Phi*sum(area(icell)));
wtL=1/(1-exp(-alffa))-1/alffa;  %weight for the trapezoid rule for time t^n
wtR=1-wtL;      %for time t^{n+1}
% % initial condition
cPrev = zeros(ncell+nedge,1);
exactCons=zeros(1,nbiter);
exactCons1=zeros(1,nbiter);
cCons=zeros(1,nbiter);
totCons=zeros(1,nbiter);
av=zeros(1,nbiter);
Lcells=cell(1,ncell);
startMMOC=0;

%% neighbors and cells around the injection well(s)
iCellNeighbors=[];
for i=1:length(pointsToTrack)
    currTris=triInvI(i,:);
    currTris(currTris==0)=[];
    for j=1:iCellCount
        if any(mainCell(currTris)==icell(j))
            iCellNeighbors=[iCellNeighbors mainCell(currTris)];
            iCellNeighbors(iCellNeighbors==icell(j))=[];
        end
    end
end
iCellNeighbors=unique(iCellNeighbors);
iCellNeighbors(iCellNeighbors==0)=[];
for m=1:nbiter
    
    %% criterion for starting MMOC
    if m>1
        checkInjNeigh=abs((cPrevPrev(iCellNeighbors)-cPrev(iCellNeighbors)));
        checkInj=abs((cPrevPrev(icell)-cPrev(icell)));
        if all(checkInjNeigh<0.0001) && all(checkInj<0.0001)  % if c is almost constant bet 2 consec time steps, use ELLAM-MMOC mix
            startMMOC=startMMOC+1;
        end
    end
    
    for i=1:ncell
        L=kappa(cg(i,1),cg(i,2))/mu(cPrev(i))*eye(2,2);
        Lcells{i}=L;%cells containing the permeability tensor L
    end
    fprintf('%i/%i\n',m,nbiter)
    %     tic;
    %%%%%%%%%%  PRESSURE/VELOCITY EQUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bp=zeros(ncell+nedge,1);
    p=zeros(ncell+nedge,1);
    nzc=0;
    for i=1:ncell
        nedge_i=size(cell_e{i},2);
        nzc = nzc + nedge_i * (1+4*nedge_i)+ncell;
    end
    IAp=zeros(nzc,1);
    JAp=zeros(nzc,1);
    VAp=zeros(nzc,1);
    %% ASSEMBLE MATRIX
    % A(IA,JA,VA) means A(IA(i),JA(i))=VA(i)
    % "pos"=position inside the vectors IA, JA, VA that store the entries of A
    posp=0;
    %% Loop over cells
    for i=1:ncell
        
        L=kappa(cg(i,1),cg(i,2))/mu(cPrev(i))*eye(2,2);
        
        
        % Compute local matrices for pressure
        Wp=local_flux_matrix(vertex(cell_v{i},:),area(i),cg(i,:),L);
        %impose zero-average condition
        for s=1:ncell
            posp=posp+1;
            IAp(posp)=i;
            JAp(posp)=s;
            VAp(posp)=area(i)*area(s);
        end
        % Loop over edges of cell
        for jj=1:size(cell_e{i},2)
            j=cell_e{i}(jj);
            
            % Inner loop over edges of cell
            for kk=1:size(cell_e{i},2)
                k=cell_e{i}(kk);
                
                % Entry cell-cell
                posp=posp+1;
                IAp(posp)=i;
                JAp(posp)=i;
                VAp(posp)=Wp(jj,kk);
                
                % Entry cell-edge k
                posp=posp+1;
                IAp(posp)=i;
                JAp(posp)=ncell+k;
                VAp(posp)=-Wp(jj,kk);
                
                % Entry edge j-cell
                posp=posp+1;
                IAp(posp)=ncell+j;
                JAp(posp)=i;
                VAp(posp)=Wp(jj,kk);
                
                % Entry edge j-edge k
                posp=posp+1;
                IAp(posp)=ncell+j;
                JAp(posp)=ncell+k;
                VAp(posp)=-Wp(jj,kk);
                
            end;
        end;
        
        % Assemble source term
        bp(i)=bp(i) + q(i);
        
    end;
    Ap=sparse(IAp(1:posp),JAp(1:posp),VAp(1:posp),ncell+nedge,ncell+nedge);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=Ap\bp;
    DarcyU=zeros(2,ncell);
    du=zeros(2,ncell);
    bc=zeros(ncell+nedge,1);
    c=zeros(ncell+nedge,1);
    nzc=0;
    posc=0;
    
    for i=1:ncell
        nedge_i=size(cell_e{i},2);
        nzc = nzc + (nedge_i+1) * (1+4*nedge_i)+1;
    end
    IAc=zeros(nzc,1);
    JAc=zeros(nzc,1);
    VAc=zeros(nzc,1);
    
    Fksigma=cell(1,ncell);
    intFlux=cell(1,ncell);
    KRDarcyU=zeros(3,sub_ncell);
    ctr=1;
    for i=1:ncell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %COMPUTATION OF GRADIENT grad p
        %FLUXES FKsigma AND
        %DARCY VELOCITY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nbe=length(cell_e{i});  %number of edges of cell i
        L=kappa(cg(i,1),cg(i,2))/mu(cPrev(i))*eye(2,2);
        [Wp,G]=local_flux_matrix(vertex(cell_v{i},:),area(i),cg(i,:),L);
        % G is the matrix with columns msigma/area * nksigma
        
        %% reconstruction via fluxes
        Fksigma{i}=Wp*(p(i)-p(ncell+cell_e{i}));
        Fksigma{i}(bdI{i})=0; %impose the strictly 0 Neumann boundary condition
        ec=edge_centers(i,cell_e,cell_v,vertex);
        DarcyU(:,i)=-(ec-ones(1,nbe)'*cg(i,:))'*Fksigma{i}/area(i);
        
        
        alfa=Phi/tstep ;
        if diffType == 1
            D=difftens(DarcyU(:,i),Phi,Phidm,Phidl,Phidt);
        else
            D=difftens_Mod(DarcyU(:,i),Phi,Phidm,Phidl,Phidt,diam(i));
        end
        KRDarcyU(:,ctr:ctr+nbe-1)=kReconstruction(Fksigma{i},area(i),vertex(cell_v{i},:),cg(i,:),vType);
        %KRDarcyU(:,i) gives the data for the reconstructed Darcy U a[x;y]+[b1;b2]
        %(a,b1,b2 respectively) on the triangular subcell i
        
        ctr=ctr+nbe;
        
        %%%%%%%%%%%%%%%% EQUATION FOR CONCENTRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% DIFFUSION MATRIX+source term
        %%
        posc=posc+1;
        IAc(posc)=i;
        JAc(posc)=i;
        if startMMOC>0
            VAc(posc)=alfa*area(i)+0.5*qp(i); %pure MMOC
        else
            VAc(posc)=alfa*area(i)+0.5*qm(i); % c for proper trap rule ,pure ELLAM
        end
        Wc=local_flux_matrix(vertex(cell_v{i},:),area(i),cg(i,:),D);
        % Loop over edges of cell
        for jj=1:size(cell_e{i},2)
            j=cell_e{i}(jj);
            
            % Inner loop over edges of cell
            for kk=1:size(cell_e{i},2)
                k=cell_e{i}(kk);
                
                % Entry cell-cell
                posc=posc+1;
                IAc(posc)=i;
                JAc(posc)=i;
                VAc(posc)=Wc(jj,kk);
                
                % Entry cell-edge k
                posc=posc+1;
                IAc(posc)=i;
                JAc(posc)=ncell+k;
                VAc(posc)=-Wc(jj,kk);
                
                % Entry edge j-cell
                posc=posc+1;
                IAc(posc)=ncell+j;
                JAc(posc)=i;
                VAc(posc)=Wc(jj,kk);
                
                % Entry edge j-edge k
                posc=posc+1;
                IAc(posc)=ncell+j;
                JAc(posc)=ncell+k;
                VAc(posc)=-Wc(jj,kk);
            end;
        end;
    end
    KRDarcyU=KRDarcyU/Phi;
    %%  Characteristic Tracking
    if startMMOC>0
        [newPointsB,newPointsF,tTrack]=complete_tracking_efficient(sub_cell_e,sub_cell_n,sub_cell_v,sub_vertex,sub_N,sub_nvert,sub_triInvVert,sIntercept,tstep,cornerPoints,pointsToTrack,pNearInj,pNearProd,triInvI,KRDarcyU);
    else
        [newPointsB,newPointsF,tTrack]=complete_tracking_efficient(sub_cell_e,sub_cell_n,sub_cell_v,sub_vertex,sub_N,sub_nvert,sub_triInvVert,sIntercept,tstep,cornerPoints,pointsToTrack,pNearInjIni,pNearProdIni,triInvI,KRDarcyU);
    end
    %% compute area of intersection of the polygons
    if startMMOC>0
        [AreaInF,AreaTf]=compAreaF(ncell,nearInjCell,subCellsIn,nxEdges,pointsToTrack, newPointsF,sub_cell_ep,pNearInj);
    else
        [AreaInB,AreaTb]=compAreaF(ncell,allCells,subCellsIn,nxEdges,pointsToTrack, newPointsB,sub_cell_ep,pNearProdIni);
        [AreaInF,AreaTf]=compAreaF(ncell,nearInjCellIni,subCellsIn,nxEdges,pointsToTrack, newPointsF,sub_cell_ep,pNearInjIni);
    end
    %AreaIn(i,j) gives the area of the intersection of cell i and tracked cell j
    invTF=cell(1,iCellCount);
    qpF=zeros(1,ncell); %ratios for the trace forward region of qp
    aTF=zeros(iCellCount,1);
    for j=1:iCellCount
        invTF{j}=find(AreaInF(:,icell(j))>0); %cells involved with the trace forward region
        invTF{j}(invTF{j}==icell(j))=[]; %remove the injection cell itself
        aTF(j)=sum(AreaInF(invTF{j},icell(j)));
        for k=1:size(invTF{j},1)
            qpF(invTF{j}(k))=AreaInF(invTF{j}(k),icell(j))/aTF(j) *qp(icell(j)); %to distribute the source evenly on region surrounding the injection well
        end
    end
    
    %%
    intCb=zeros(1,ncell);
    intCf=zeros(1,ncell);
    for i=1:ncell
        intCb(i)=AreaInB(:,i)'*cPrev(1:ncell); %integral of c, backtracked
        intCf(i)=AreaInF(i,:)*cPrev(1:ncell); %integral of c, forward tracked
        if startMMOC>0 %MMOC
            bc(i)=bc(i) + intCf(i)*Phi/tstep-0.5*cPrev(i)*qp(i)+qp(i);
        else % ELLAM
            bc(i)=bc(i)-0.5*qm(i)*cPrev(i) + intCb(i)*Phi/tstep+wtR*qp(i)+wtL*(qpF(i))*(1-exp(-alffa))+wtL*exp(-alffa)*qp(i); %c for proper trap rule
        end
    end
    
    Ac=sparse(IAc(1:posc),JAc(1:posc),VAc(1:posc),ncell+nedge,ncell+nedge);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c=Ac\bc;
    
    %% measures how well c satisfies global mass conservation
    for i=1:ncell
        av(m)=av(m)+cPrev(i)*area(i);
        exactCons(m)=exactCons(m)+tstep*qp(i)-wtL*tstep*qm(i)*cPrev(i)-wtR*tstep*qm(i)*c(i);
        cCons(m)=cCons(m)+Phi*c(i)*area(i)-Phi*cPrev(i)*area(i)-tstep*qp(i)+wtL*tstep*qm(i)*cPrev(i)+wtR*tstep*qm(i)*c(i);
    end
    totCons(m)=sum(cCons)/sum(exactCons);
    exactCons1(m)=sum(exactCons);
    av(m)=av(m)/sum(area);
    
    %% re-initialize c
    cPrevPrev=cPrev;
    cPrev=c;
end
eCons=cCons./exactCons1;
toc
% save('c_Hexa.mat','c','av','pPerEdge','castest','mesh');
% write_solution_vtk(c,'solution_concentration',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg,area);
% write_solution_vtk(p,'solution_pressure',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
% write_solution_vtk(zeros(size(c)),'grid',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
