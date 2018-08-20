function D = interiorkjli(x,a,uk,ul,c,fi,fj)
%Calculate donation values given interior solution given kj and li are interior.


D(1) = (2*pi)^(-0.5)*a*(x(1))^(a-1)*exp(-0.5*((x(1))^a+fi^a-(x(2))^a-fj^a)^2)*ul-c;
D(2) = -(2*pi)^(-0.5)*a*(x(2))^(a-1)*exp(-0.5*((x(1))^a+fi^a-(x(2))^a-fj^a)^2)*uk-c;

end

