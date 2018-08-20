function D = interiorlj(x,a,ul,kbar,c,fi,fj)
%Calculate donation values given ki=kbar and lj is interior.


D = -(2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((kbar)^a+fi^a-(x)^a-fj^a)^2)*ul-c;

end

