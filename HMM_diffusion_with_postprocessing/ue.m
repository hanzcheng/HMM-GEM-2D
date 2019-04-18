function [uv,avu]=ue(z)
global testcase
x=z(:,1);y=z(:,2);
if testcase == 1 
uv=x;
avu=0.5;
elseif testcase == 2
uv=x .^ 2 + y .^ 2;
avu=2/3;
elseif testcase == 3
uv=(x-x.^2).*(y-y.^2);
avu=1/36;
elseif testcase == 4
uv=sin(pi*x) .* sin(pi*y);
avu=4/(pi)^2;
elseif testcase == 5
uv=(x - exp(20000*(x - 1))) .* (y.^2 - exp(30000*(y - 1)));
avu=1;
elseif testcase == 6
uv=cos(pi*x) .* cos(pi*y);
avu=0;
end
%uv=x*sin(pi*x)*y*cos((pi*y)/2);
