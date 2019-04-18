function [ux,uy,H]=derivu(z)
global testcase
x=z(1);
y=z(2);
H=zeros(2,2);

if testcase == 1
 ux=1;
 uy=0;
elseif testcase == 2
    ux=2*x;
    uy=2*y;
    H(1,1)=2; %u_xx
    H(1,2)=0; %u_xy
    H(2,2)=2; %u_yy
    H(2,1)=H(1,2); %u_yx
elseif testcase == 3
    ux=(1-2*x)*y*(1-y);
    uy=x*(1-x)*(1-2*y);
    H(1,1)=-2*y*(1-y);
    H(1,2)=(1-2*x)*(1-2*y);
    H(2,2)=-2*x*(1-x);
    H(2,1)=H(1,2);
elseif testcase ==4
    ux= pi*cos(pi*x)*sin(pi*y);
    uy=pi*cos(pi*y)*sin(pi*x);
    H(1,1)=-pi^2*sin(pi*x)*sin(pi*y);
    H(1,2)=pi^2*cos(pi*x)*cos(pi*y);
    H(2,1)=H(1,2);
    H(2,2)=-pi^2*sin(pi*y)*sin(pi*x);
    elseif testcase ==6
    ux= -pi*sin(pi*x)*cos(pi*y);
    uy= -pi*sin(pi*y)*cos(pi*x);
    H(1,1)=-pi^2*cos(pi*x)*cos(pi*y);
    H(1,2)=pi^2*sin(pi*x)*sin(pi*y);
    H(2,1)=H(1,2);
    H(2,2)=-pi^2*cos(pi*y)*cos(pi*x);

end

