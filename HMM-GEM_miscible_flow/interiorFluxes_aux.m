function intFlux = interiorFluxes_aux(Fksigma,area,vertex,center)
solMethod = "mat"; % compute the q_i via matrix inversion in matlab
solMethod = "exp"; %other option: "exp" for explicit expression
s=size(Fksigma,1); %gives number of edges since Fksigma is edge fluxes
P=eye(s,s);
b=zeros(s,1);
aT=zeros(s,1);
intFlux = zeros(s,1); %interior fluxes
oN = zeros(2,s); %outward normal along each side
v1=vertex(1,:);
v2=vertex(2,:);
v3=center;
aT(1)=0.5*det([v1-v3;v2-v3]);
F_k = zeros(1,2);
F_star = zeros(s,1);
beta_i = zeros(s,1);
q = zeros(s,1);
diam = zeros(s,1); % diameter of subcells
edge_len = zeros(s,1); %edge lengths
% construct a cell-wise consistent flux
for i=1:s
    v1 = vertex(i,:);
    v2 = vertex(i+1,:);
    v3 = center;
    edge_len(i) = norm(v3-v2);
    diam(i) = max(norm(v2-v1),norm(v3-v1));
    diam(i) = max(diam(i), norm(v3-v2));
    oN(:,i) = [0 1; -1 0] * (v3-v2)'/norm(v3-v2);
    vMid = 0.5*(v1+v2);
    aT(i)=0.5*det([v1-v3;v2-v3]); %Area of Triangle i
    F_k = F_k + Fksigma(i)* (vMid-v3);
    if i<s
        P(i,i+1) = -1;
    else
        P(i,1) = -1;
    end
end
F_k = F_k / area;
for i=1:s
    F_star(i) = edge_len(i)*F_k*oN(:,i);
end
cons_ratio = zeros(s,1);
for i=1:s
    if i<s
        beta_i(i) = 2*edge_len(i)/(diam(i)+diam(i+1));
    else
        beta_i(i)= 2*edge_len(i)/(diam(i)+diam(1));
    end
    %% Ratio for C-velocity
    if i<s
        cons_ratio(i) = (aT(i)+aT(i+1))/(2*area);
    else
        cons_ratio(i) = (aT(i)+aT(1))/(2*area);
    end
    %% to make this equal to C-velocity, uncomment latter statement
    %     beta_i(i) = 1/cons_ratio(i);
    if i>1
        b(i)=-F_star(i)-Fksigma(i)+aT(i)/area*sum(Fksigma)+F_star(i-1);
    else
        b(i)=-F_star(i)-Fksigma(i)+aT(i)/area*sum(Fksigma)+F_star(s);
    end
end
%% explicit expression for q_i
if solMethod=="exp"
    num_qs = 0;
    den_qs = sum(1./beta_i);
    for i=1:s-1
        num_qs = num_qs + sum(b(i)./beta_i(1:s-1));
    end
    q(s) = -1/beta_i(s)*num_qs/den_qs;
    q(1) = (beta_i(s)*q(s)+b(1))/beta_i(1);
    for i=2:s-1
        q(i) = (q(i-1)*beta_i(i-1)+b(i))/beta_i(i);
    end
    %% solve for q_i using matlab's matrix inversion
else
    D = diag(beta_i);
    A = P'*D;
    A(s,:) = 1;
    b(s) = 0;
    q = A\b;
end
intFlux = F_star + beta_i .* q;
end