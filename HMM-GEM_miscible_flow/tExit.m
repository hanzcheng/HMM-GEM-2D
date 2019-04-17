function [t_Exit,nextTri]=tExit(currTri,sub_cell_e,sub_cell_n,sIntercept,KRDarcyU,t_left,pointsToTrack,tracking)
xToTrack=pointsToTrack(:,1);
yToTrack=pointsToTrack(:,2);
nbPts=size(pointsToTrack,1);
t_Exit=zeros(nbPts,1);
nextTri=zeros(nbPts,1);
tol=1e-8;
for i=1:nbPts
    timExit=Inf;
    x=xToTrack(i);
    y=yToTrack(i);
    triNow=currTri(i); %triangles involved for the ith point
    a=KRDarcyU(1,triNow);
    b1=KRDarcyU(2,triNow);
    b2=KRDarcyU(3,triNow);
    m=sIntercept(1,sub_cell_e{triNow});
    d=sIntercept(2,sub_cell_e{triNow});
    onEdge=m*x+d-y;
    onEdge(m==Inf)=x-d(m==Inf);
    notOnEdge=find(abs(onEdge)>tol);
    temptime=zeros(size(notOnEdge));
    for k=1:size(notOnEdge,2)
        if a==0
            if tracking==0 %back track
                if(abs(m(notOnEdge(k)))==Inf)
                    temptime(k)=t_left(i)+(d(notOnEdge(k))-x)/b1;
                else
                    temptime(k)=t_left(i)+(m(notOnEdge(k))*x+d(notOnEdge(k))-y)/(b2-m(notOnEdge(k))*b1);
                end
                temptime(k)=t_left(i)-temptime(k);
            elseif tracking==1 %forward track
                if(abs(m(notOnEdge(k)))==Inf)
                    temptime(k)=(d(notOnEdge(k))-x)/b1;
                else
                    temptime(k)=(m(notOnEdge(k))*x+d(notOnEdge(k))-y)/(b2-m(notOnEdge(k))*b1);
                end
            end
            
        else
            if tracking==0
                if(abs(m(notOnEdge(k)))==Inf)
                    temptime(k)=t_left(i)+ (1/a) * log((b1+a*d(notOnEdge(k)))/(a*x+b1));
                else
                    temptime(k)=t_left(i)+ (1/a) * log((b2-m(notOnEdge(k))*b1+d(notOnEdge(k))*a)/(a*y+b2-m(notOnEdge(k))*(a*x+b1)));
                end
                temptime(k)=t_left(i)-temptime(k);
            elseif tracking==1
                if(abs(m(k))==Inf)
                    temptime(k)=(1/a) * log((b1+a*d(notOnEdge(k)))/(a*x+b1));
                else
                    temptime(k)=(1/a) * log((b2-m(notOnEdge(k))*b1+d(notOnEdge(k))*a)/(a*y+b2-m(notOnEdge(k))*(a*x+b1)));
                end
            end
        end
        if temptime(k)<tol || ~isreal(temptime(k))
            temptime(k)=Inf;
        end
        if abs(timExit)>abs(temptime(k))
            timExit=temptime(k);
            eExit=notOnEdge(k);
            nextTri(i)=sub_cell_n{currTri(i)}(eExit);
        end
    end
    t_Exit(i)=timExit;
end
t_Exit(t_Exit>t_left)=t_left(t_Exit>t_left);
nextTri(nextTri==0)=currTri(nextTri==0);
