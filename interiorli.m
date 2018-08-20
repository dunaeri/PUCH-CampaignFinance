function D = interiorli(x,a,ul,kbar,c,fi,fj)
%Calculate donation values given kj=kbar and li is interior.


D = (2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((x)^a+fi^a-(kbar)^a-fj^a)^2)*ul-c;

end

