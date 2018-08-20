function D = interiorki2(x,a,uk,c,fi,fj)
%Calculate donation values given ki is interior and lj=0.


D = (2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((x)^a+fi^a-fj^a)^2)*uk-c;

end

