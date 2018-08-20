function D = interiorlj2(x,a,ul,c,fi,fj)
%Calculate donation values given ki=0 and lj is interior.


D = -(2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*(fi^a-(x)^a-fj^a)^2)*ul-c;

end

