function D = interiorli2(x,a,ul,c,fi,fj)
%Calculate donation values given kj=0 and li is interior.


D = (2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((x)^a+fi^a-fj^a)^2)*ul-c;

end