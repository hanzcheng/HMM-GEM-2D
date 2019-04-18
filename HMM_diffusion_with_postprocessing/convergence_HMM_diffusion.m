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

%% load a sequence of meshes to measure convergence
%meshes={'mesh2_1.mat';'mesh2_2.mat';'mesh2_3.mat';'mesh2_4.mat'};%'mesh2_5.mat'};
meshes={'mesh2_1.mat';'mesh2_2.mat';'mesh2_3.mat'};
%meshes={'mesh3_1.mat';'mesh3_2.mat';'mesh3_3.mat';'mesh3_4.mat'};
% meshes={'mesh3_1.mat';'mesh3_2.mat';'mesh3_3.mat'};
% meshes={'mesh4_1.mat'};
% meshes={'hexa1_2.mat'};
nbmeshes=size(meshes,1);

for mm=1:nbmeshes
    tic;
    % Load mesh
    loadmesh=strcat('load ..//..//matlab_meshes/',meshes{mm});
    disp(loadmesh);
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
    
    for i=1:ncell
        nedge_i=size(cell_e{i},2);
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
        
        % Assemble the source term f1
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
    %% Computation mesh size and max error
    h(mm) = max(abs(diam));
    if bccase == 1 %pure Neumann
        [uex,avu]=ue(center);
        
        Linfe(mm)=max(abs((u(1:ncell)-ue(cg)+avu)));
        pLinfe(mm)=Linfe(mm)/max(ue(cg)-avu);
        L2e(mm)=sqrt(sum(area.*(u(1:ncell)-ue(cg)+avu).^2));
        pL2e(mm)=L2e(mm)/sqrt(sum(area.*(ue(cg)-avu).^2));
    else
        Linfe(mm)=max(abs((u(1:ncell)-ue(cg))));
        pLinfe(mm)=Linfe(mm)/max(ue(cg));
        L2e(mm)=sqrt(sum(area.*(u(1:ncell)-ue(cg)).^2));
        pL2e(mm)=L2e(mm)/sqrt(sum(area.*(ue(cg)).^2));
    end
    
    %% compute the flux and flux error
    Fksigma=cell(1,ncell);
    F_exact = cell(1,ncell);
    maxErr=0; % maximum error in flux approximation
    for i=1:ncell
        nbe=length(cell_e{i});
        [Wp,G]=local_flux_matrix(vertex(cell_v{i},:),area(i),center(i,:),cg(i,:));
        Fksigma{i}=Wp*(u(i)-u(ncell+cell_e{i}));
        edge_centres = 0.5*(vertex(cell_v{i}(1:nbe),:)+vertex(cell_v{i}(2:nbe+1),:));
        edge_lengths = (vertex(cell_v{i}(1:nbe),:)-vertex(cell_v{i}(2:nbe+1),:)).^2;
        edge_lengths = sqrt(edge_lengths(:,1)+edge_lengths(:,2));
        F_exact{i} = zeros(nbe,1);
        for j=1:nbe
            v1=vertex(cell_v{i}(j),:);
            v2=vertex(cell_v{i}(j+1),:);
            N = (v2-v1) * [0 -1;1 0];
            msigma=sqrt(sum(N.^2,2));
            % Normalisation of N
            
            N=N./[msigma msigma];
            [ux, uy]= derivu(edge_centres(j,:));
            gradu = [ux;uy];
            L=lambda(edge_centres(j,:));
            F_exact{i}(j)=-N*L*gradu;
            Fksigma{i}(j) = Fksigma{i}(j)/msigma;
            maxErr=max(maxErr,abs(F_exact{i}(j)-Fksigma{i}(j)));
        end
    end
end
pLinfe
pL2e
h

%% convergence in L2 and Linf
for ml=1:nbmeshes-1 
    oclinf(ml)=log(Linfe(ml)/Linfe(ml+1))/log(h(ml)/h(ml+1));
    ocl2(ml)=log(L2e(ml)/L2e(ml+1))/log(h(ml)/h(ml+1));
end

if (nbmeshes>1)
    oclinf
    ocl2
end
% also measure error in the gradient
[L2u,L2gradu]=compute_errors(u,area,center,vertex,cell_v,cell_e);
write_solution_vtk(u,'solution',ncell,nedge,nvert,cell_v,vertex)
%[time_loadmesh time_assembly time_solution]

