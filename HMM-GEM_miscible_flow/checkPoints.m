function posSame=checkPoints(pointsToCheck,sub_vertex)
%posSame(i) = 0 means that pointsToCheck(i) does not coincide with any sub vertex
%posSame(i) neq 0 means that it coincides with the posSame(i)th sub vertex
nbPtsCheck=size(pointsToCheck,1);
nbPtsMain=size(sub_vertex,1);
posSame=zeros(nbPtsCheck,1);
thres=1e-5;
for i=1:nbPtsCheck
    xydiff=[pointsToCheck(i,1)*ones(nbPtsMain,1) pointsToCheck(i,2)*ones(nbPtsMain,1)]-sub_vertex;
    normXy= sqrt(sum(abs(xydiff).^2,2));
    [~,posMin]=min(normXy);
    if norm(sub_vertex(posMin,:))~=0
        if norm(xydiff(posMin,:))/norm(sub_vertex(posMin,:))<=thres
            [~,posSame(i)]=min(normXy); %gives the location of the sub_vertex nearest to pointsToCheck(i,:)
        end
    else
        if min(normXy)<=thres
        [~,posSame(i)]=min(normXy); %gives the location of the sub_vertex nearest to pointsToCheck(i,:)
        end
    end
end
