function D = interiorki(x,a,uk,lbar,c,fi,fj)
%Calculate donation values given ki is interior and lj=lbar.


D = (2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((x)^a+fi^a-(lbar)^a-fj^a)^2)*uk-c;

end

