function D = normal(x,a,uk,ul)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


D(1) = (2*pi)^(-0.5)*a*x(1)^(a-1)*exp(-0.5*(x(1)^a-x(2)^a)^2)*uk-1;
D(2) = -(2*pi)^(-0.5)*a*x(2)^(a-1)*exp(-0.5*(x(1)^a-x(2)^a)^2)*ul-1;

end

