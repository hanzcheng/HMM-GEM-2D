function currTri=detTri(TriInv,prevTri,sub_cell_n,sub_cell_v,sub_vertex,sub_N,KRDarcyU,pointsToTrack,track,locCorner)
xToTrack=pointsToTrack(:,1);
yToTrack=pointsToTrack(:,2);
nbPts=size(pointsToTrack,1);
currTri=zeros(nbPts,1);
tol=1e-10;
for i=1:nbPts

    x=xToTrack(i);
    y=yToTrack(i);
    posTriNow=unique(TriInv(i,:));
    triNow=posTriNow(posTriNow~=0); %triangles involved for the ith point
    nTriNow=length(triNow); %to know how many triangles are involved
    n1=[];
    n2=[];
    n3=[];
    if nTriNow==1
        currTri(i)=triNow;
    elseif nTriNow==2  && ~any(locCorner==i)%this means that we are looking at a point that lies on an edge not on the corners
        T1=triNow(1);
        Tp=prevTri(i);
        
        a=KRDarcyU(1,T1);
        b1=KRDarcyU(2,T1);
        b2=KRDarcyU(3,T1);
        locDarcyU=a*[x y]+[b1 b2];
        if any(triNow==sub_cell_n{T1}(1))
            n1=1;
        end
        if any(triNow==sub_cell_n{T1}(2))
            n2=2;
        end
        if any(triNow==sub_cell_n{T1}(3))
            n3=3;
        end
        Linv=[n1 n2 n3]';
        oN=sub_N{T1}(Linv,:);%to get the outer normal of the line on which the point lies on
        uDotN=oN*locDarcyU';
        if T1==Tp
            currTri(i)=triNow(2);
        elseif uDotN>0
            if track==0
                currTri(i)=T1;
            elseif track==1
                currTri(i)=triNow(2);
            end
        elseif uDotN<0
            if track==0
                currTri(i)=triNow(2);
            elseif track==1
                currTri(i)=T1;
            end
        else
            currTri(i)=T1;
        end
    else %we are looking at a point that lies on a vertex.
        
        ctr=1;
        s1=zeros(nTriNow,1);
        s2=zeros(nTriNow,1);
        dfZ=zeros(nTriNow,1);
        Tp=prevTri(i); %previous triangle
        
        if Tp~=0
            a=KRDarcyU(1,Tp);
            b1=KRDarcyU(2,Tp);
            b2=KRDarcyU(3,Tp);
        end
        
        while ctr<=nTriNow && currTri(i)==0
            Linv=[];
            T1=triNow(ctr);
            if Tp==0
                a=KRDarcyU(1,T1);
                b1=KRDarcyU(2,T1);
                b2=KRDarcyU(3,T1);
            end
            locDarcyU=a*[x y]+[b1 b2];
            if track==0
                locDarcyU=-locDarcyU;
            end
            locDarcyU(abs(locDarcyU)<tol)=0;
            nM1=norm([x y]-sub_vertex(sub_cell_v{T1}(1),:));
            nM2=norm([x y]-sub_vertex(sub_cell_v{T1}(2),:));
            nM3=norm([x y]-sub_vertex(sub_cell_v{T1}(3),:));
            if nM1<nM2
                if nM1<nM3
                    posV=1; %position of the vertex
                else
                    posV=3;
                end
            else
                if nM2<nM3
                    posV=2;
                else
                    posV=3;
                end
            end
            if posV==1
                Linv=[1;3];
            elseif posV==2
                Linv=[1;2];
            elseif posV==3
                Linv=[2;3];
            end
            oN=sub_N{T1}(Linv,:);
            v=oN*[0 1; -1 0];
            v1=v(1,:);
            v2=v(2,:);
            if posV==1
                v2=-v2;
            elseif posV==2
                v1=-v1;
            elseif posV==3
                v1=-v1;
            end
            
            if T1 ~= Tp
                s1(ctr)=-det([locDarcyU' v1'])*det([v1' v2']);
                s2(ctr)=det([locDarcyU' v2'])*det([v1' v2']);
                sg1=sign(s1(ctr));
                sg2=sign(s2(ctr));
                if sg1>=0 && sg2>=0
                    currTri(i)=T1;
                end
                if abs(s1(ctr))<tol
                    sg1=0;
                end
                if abs(s2(ctr))<tol
                    sg2=0;
                end
                if sg1>=0 && sg2>=0
                    currTri(i)=T1;
                end
                if sg1<0
                    dfZ(ctr)=dfZ(ctr)+abs(s1(ctr));
                end
                if sg2<0
                    dfZ(ctr)=dfZ(ctr)+abs(s2(ctr));
                end
            end
            ctr=ctr+1;
        end
        if nTriNow==2 && currTri(i)==0 %corner point
            T1=triNow(1);
            T2=triNow(2);
            if Tp==T1
                currTri(i)=T2;
            elseif Tp==T2
                currTri(i)=T1;
            end
        end
        dfZ(triNow==Tp)=Inf;
        if currTri(i)==0
            [~,index]=min(dfZ);
            currTri(i)=triNow(index);
        end
    end
end
end