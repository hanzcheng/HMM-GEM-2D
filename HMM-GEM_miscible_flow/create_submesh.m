function [sub_ncell,sub_nvert,sub_nedge,sub_vertex,sub_cell_v,sub_cell_n,sub_cell_e,sub_N,sub_triInvVert,mainCell,subCellsIn]=create_submesh(ncell,nvert,vertex,cg,cell_v,cell_e)
sub_vertex=[vertex(:,:);cg(:,:)]; %vertices of the submesh
sub_ncell=0; %number of cells in the submes
subCellsIn=cell(1,ncell); % subCellsIn{i} gives the numbers of the sub cells (triangles) in cell i
for i=1:ncell
    nbe=length(cell_e{i});
    sub_ncell=sub_ncell+nbe;
end
sub_cell_v=cell(1,sub_ncell); %numbers of the vertices of the submesh
sub_cell_n=cell(1,sub_ncell); %numbers of the neighboring cells of the submesh
sub_cell_e=cell(1,sub_ncell); %numbers of the edges of the submesh
sub_N=cell(1,sub_ncell); %gives the unit outer normal to the edges of each subcell
 %to access the outward normal along edge j of subcell i, we use
 %sub_N{i}(j,:)
ctr=1;
sub_nvert=length(sub_vertex);
sub_triInvVert=zeros(sub_nvert,sub_nvert); 
%triInvVert(i,:) gives triangles involved with vertex i
mainCell=zeros(1,sub_ncell);
%mainCell(i) gives the main cell where triangle i came from
edgeNum=1; %number of the edges
Edge_Connections=zeros(sub_nvert,sub_nvert);
Edge_Num=zeros(3,sub_ncell);
for i=1:ncell
    nbe=length(cell_e{i});
    subCellsIn{i}=zeros(1,nbe);
    for j=1:nbe
        subCellsIn{i}(j)=ctr;
        mainCell(ctr)=i;
        sub_cell_v{ctr}=[cell_v{i}(j), cell_v{i}(j+1), nvert+i, cell_v{i}(j)];
        z1=find(sub_triInvVert(cell_v{i}(j),:)==0);
        sub_triInvVert(cell_v{i}(j),z1(1))=ctr;
        z2=find(sub_triInvVert(cell_v{i}(j+1),:)==0);
        sub_triInvVert(cell_v{i}(j+1),z2(1))=ctr;
        z3=find(sub_triInvVert(nvert+i,:)==0);
        sub_triInvVert(nvert+i,z3(1))=ctr;
        if Edge_Connections(cell_v{i}(j),cell_v{i}(j+1))==0
            Edge_Connections(cell_v{i}(j),cell_v{i}(j+1))=edgeNum;
            Edge_Connections(cell_v{i}(j+1),cell_v{i}(j))=edgeNum;
            edgeNum=edgeNum+1;
        end
        if Edge_Connections(cell_v{i}(j+1),nvert+i)==0
            Edge_Connections(cell_v{i}(j+1),nvert+i)=edgeNum;
            Edge_Connections(nvert+i,cell_v{i}(j+1))=edgeNum;
            edgeNum=edgeNum+1;
        end
        if Edge_Connections(cell_v{i}(j),nvert+i)==0
            Edge_Connections(cell_v{i}(j),nvert+i)=edgeNum;
            Edge_Connections(nvert+i,cell_v{i}(j))=edgeNum;
            edgeNum=edgeNum+1;
        end
        sub_cell_e{ctr}=[Edge_Connections(cell_v{i}(j),cell_v{i}(j+1)), Edge_Connections(cell_v{i}(j+1),nvert+i),  Edge_Connections(cell_v{i}(j),nvert+i)];
        Edge_Num(:,ctr)=[Edge_Connections(cell_v{i}(j),cell_v{i}(j+1)); Edge_Connections(cell_v{i}(j+1),nvert+i);  Edge_Connections(cell_v{i}(j),nvert+i)];
        sub_cell_n{ctr}=zeros(1,3);
        v1=sub_vertex(sub_cell_v{ctr}(1:3),:);
        v2=sub_vertex(sub_cell_v{ctr}(2:4),:);
        %to get the outward normal
        N = (v2-v1) * [0 -1;1 0];
        % Computation of the length of the edge
        msigma=sqrt(sum(N.^2,2));
        % Normalisation of N
        N=N./[msigma msigma];
        sub_N{ctr}=N;
        ctr=ctr+1;
        
    end
    
end
sub_nedge=edgeNum-1;
for i=1:sub_ncell
    for j=1:3
        [~,q]=find(Edge_Num==sub_cell_e{i}(j));
        q(q==i)=[];
        if ~isempty(q)
            sub_cell_n{i}(j)=q;
        end
    end
end