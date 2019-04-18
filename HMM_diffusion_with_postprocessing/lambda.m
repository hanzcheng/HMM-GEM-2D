function [L,dL]=lambda(z) %diffusion matrix
x=z(:,1);
y=z(:,2);
global diffmat
L=eye(2,2);
dL=zeros(2,2);
if diffmat ==1
    L=eye(2,2);
    dL=zeros(2,2);
elseif diffmat == 2
    L(1,1)=1.5;
    L(2,2)=1.5;
    L(1,2)=0.5;
    dL(1,1)=0; %1st row of dL gives partial derivatives wrt x
    dL(1,2)=0;
    dL(2,1)=0; %2nd row of dL gives partial derivatives wrt y
    dL(2,2)=0;
elseif diffmat == 3
    L(1,1)=1;
    L(2,2)=1000;
    L(1,2)=0;
    dL(1,1)=0; %1st row of dL gives partial derivatives wrt x
    dL(1,2)=0;
    dL(2,1)=0; %2nd row of dL gives partial derivatives wrt y
    dL(2,2)=0;
elseif diffmat == 4
    L(1,1)=1;
    L(2,2)=10^-1;
    R=[cosd(40),-cosd(40);sind(40),sind(40)];
    L=R*L/R;
    dL(1,1)=0; %1st row of dL gives partial derivatives wrt x
    dL(1,2)=0;
    dL(2,1)=0; %2nd row of dL gives partial derivatives wrt y
    dL(2,2)=0;
elseif diffmat == 5
    nu=10^-3;
    adj=0.1;
    x=x+adj;
    y=y+adj;
    L(1,1)=nu* (x^ 2)+y^ 2;
    L(1,2)=(nu - 1)* x * y;
    L(2,1)=L(1,2);
    L(2,2)=nu* (y^2)+x^2;
    L=1 / (x^2 + y^2) *L;
    dL(1,1)=2*x*y^2*(nu -1)/(x^2+y^2);
    dL(1,2)=(nu-1)*y*(y^2-x^2)/(x^2+y^2) ;
    dL(2,1)=(nu-1)*x*(x^2-y^2)/(x^2+y^2)  ;
    dL(2,2)=2*x^2*y*(nu -1)/(x^2+y^2);
    dL=1/(x^2 + y^2) *dL;
end

L(2,1)=L(1,2);
