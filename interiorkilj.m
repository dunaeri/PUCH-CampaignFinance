function D = interiorkilj(x,a,uk,ul,c,fi,fj)
%Calculate donation values given interior solution given ki and lj are interior.


D(1) = (2*pi)^(-0.5)*a*(x(1))^(a-1)*exp(-0.5*((x(1))^a+fi^a-(x(2))^a-fj^a)^2)*uk-c;
D(2) = -(2*pi)^(-0.5)*a*(x(2))^(a-1)*exp(-0.5*((x(1))^a+fi^a-(x(2))^a-fj^a)^2)*ul-c;

end

