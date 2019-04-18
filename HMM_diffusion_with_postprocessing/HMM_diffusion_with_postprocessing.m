%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HMM for  - div (lambda grad u)   = f1   in  \Omega
%          Dirichlet or Neumann  conditions   on  \partial\Omega_D
% Here first numbering edges and then cells.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long
global testcase
testcase = 6;  % different test cases correspond to different exact sol. 'u'
global bccase
bccase = 1; %boundary condition 0 refers to pure dirichlet conditions while
% boundary condition 1 refers to pure Neumann conditions
% boundary condition 2: mixed Dirichlet-Neumann (Neumann at x=1, Dirichlet
% elswhere)
global diffmat
diffmat = 1; %diffmat 1 means that lambda is the identity matrix
%diffmat 2 onwards corresponds to other diffusion tensors

%mesh
meshes={'mesh2_3.mat'};
% meshes={'mesh4_1.mat'};
% meshes={'hexa1_2.mat'};

    tic;
    % Load mesh
    loadmesh=strcat('load ..//..//matlab_meshes/',meshes{1});
    eval(loadmesh);
    cg=gravity_centers(ncell,cell_v,vertex,area);
    disp('mesh loaded');
    time_loadmesh=toc;
    
    tic;
    % Matrix, unknown and RHS
    %		The unknows are the cells (1:ncell) and the edges (ncell+1:ncell+nedge)
    b=zeros(ncell+nedge,1);
    u=zeros(ncell+nedge,1);
    nz=0;
    
    for i=1:ncell
        cell_n{i}=SetBC(cell_n{i},vertex(cell_v{i},:));
    end

    totEdges = 0; %total number of edges
    for i=1:ncell
        nedge_i=size(cell_e{i},2); % number of edges in cell i
        totEdges = totEdges + nedge_i;
        nz = nz + nedge_i * (1+4*nedge_i);
    end
    if bccase==1
        nz = nz+ncell*ncell;
    end
    IA=zeros(nz,1);
    JA=zeros(nz,1);
    VA=zeros(nz,1);
    %% ASSEMBLE MATRIX
    % A(IA,JA,VA) means A(IA(i),JA(i))=VA(i)
    % "pos"=position inside the vectors IA, JA, VA that store the entries of A
    pos=0;
    %% Loop over cells
    for i=1:ncell
        
        % Compute local matrices
        W=local_flux_matrix(vertex(cell_v{i},:),area(i),center(i,:),cg(i,:));
        % the inputs come from the generated mesh
        if bccase ==1
            for s=1:ncell
                pos=pos+1;
                IA(pos)=i;
                JA(pos)=s;
                VA(pos)=area(i)*area(s);
            end
        end
        % Loop over edges of cell
        for jj=1:size(cell_e{i},2)
            j=cell_e{i}(jj);
            % Dirichlet boundary condition
            if (cell_n{i}(jj)==0)
                pos=pos+1;
                IA(pos)=ncell+j;
                JA(pos)=ncell+j;
                VA(pos)=1;
            end
            
            % Inner loop over edges of cell
            for kk=1:size(cell_e{i},2)
                k=cell_e{i}(kk);
                
                % Entry cell-cell
                pos=pos+1;
                IA(pos)=i;
                JA(pos)=i;
                VA(pos)=W(jj,kk);
                
                % Entry cell-edge k
                pos=pos+1;
                IA(pos)=i;
                JA(pos)=ncell+k;
                VA(pos)=-W(jj,kk);
                
                % If j is not a Dirichlet boundary edge
                if (cell_n{i}(jj)~=0)
                    % Entry edge j-cell
                    pos=pos+1;
                    IA(pos)=ncell+j;
                    JA(pos)=i;
                    VA(pos)=W(jj,kk);
                    
                    % Entry edge j-edge k
                    pos=pos+1;
                    IA(pos)=ncell+j;
                    JA(pos)=ncell+k;
                    VA(pos)=-W(jj,kk);
                end
            end;
        end;
        
        % Assemble source term
        b(i)=b(i) + area(i)*source_term(cg(i,:));
        Dbe=find(cell_n{i}==0); % to determine which edge is a Dirichlet boundary edge
        Nbe=find(cell_n{i}==-1); % to determine which edge is a Neumann boundary edge
        Dvb1=vertex(cell_v{i}(Dbe),:); %1st vertex of the boundary edge(s)
        Dvb2=vertex(cell_v{i}(Dbe+1),:); %2nd vertex of the boundary edge(s)
        Nvb1=vertex(cell_v{i}(Nbe),:); %1st vertex of the boundary edge(s)
        Nvb2=vertex(cell_v{i}(Nbe+1),:); %2nd vertex of the boundary edge(s)
        cDvb=(Dvb1+Dvb2)/2;
        
        %% Dirichlet Boundary Conditions
        if Dbe~=0
            dbc=ue(cDvb);
            b(ncell+cell_e{i}(Dbe))=dbc;
        end
        %% Neumann Boundary Conditions
        %non-homogeneous
        if Nbe~=0
            L=lambda(cg(i,:));
            g=Neumann_BC(Nvb1,Nvb2,L);
            b(ncell+cell_e{i}(Nbe))=b(ncell+cell_e{i}(Nbe))-g;
        end
    end;
    
    A=sparse(IA(1:pos),JA(1:pos),VA(1:pos),ncell+nedge,ncell+nedge);
    
    disp('Matrix assembled')
    time_assembly=toc;
    
    tic;
    % Solution
    u=A\b;
    disp('Solution computed')
    time_solution=toc;
    
    %% compute the fluxes and reconstruct an RT0 velocity field
    
    Fksigma=cell(1,ncell); %fluxes
    
    KRDarcyU=zeros(3,totEdges); % RT0 velocity field
    vType = "A"; %choice of reconstruction
%     vType = "C";
%     vType = "KR";
    ctr = 1; % keeps track of which sub-cells the velocity field is being constructed on
    for i=1:ncell
        nbe=length(cell_e{i});
        [Wp,G]=local_flux_matrix(vertex(cell_v{i},:),area(i),center(i,:),cg(i,:));
        Fksigma{i}=Wp*(u(i)-u(ncell+cell_e{i}));
        KRDarcyU(:,ctr:ctr+nbe-1)=kReconstruction(Fksigma{i},area(i),vertex(cell_v{i},:),cg(i,:),vType);
        ctr = ctr+nbe;
    end

    % KRDarcyU is an RT0 velocity defined piecewise on the following submesh:
    [sub_ncell,sub_nvert,sub_nedge,sub_vertex,sub_cell_v,sub_cell_n,sub_cell_e,sub_N,sub_triInvVert,mainCell,subCellsIn]=create_submesh(ncell,nvert,vertex,cg,cell_v,cell_e);
    % in particular, KRDarcyU(:,i) is the RT0 velocity for the i'th
    % cell (triangle) defined in this mesh 
    % This RT velocity is useful, for example, when performing a
    % characteristic tracking, as implemented in the miscible flow model
write_solution_vtk(u,'solution',ncell,nedge,nvert,cell_v,vertex)
%[time_loadmesh time_assembly time_solution]

