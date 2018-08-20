function D = interiordj(x,a,uk,c,fi,fj)
%Calculate donation values given both donors support candidate j.


D = -(2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*(fi^a-(x)^a-fj^a)^2)*uk-c;

end

