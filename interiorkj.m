function D = interiorkj(x,a,uk,lbar,c,fi,fj)
%Calculate donation values given kj=interior and li=lbar.


D = -(2*pi)^(-0.5)*a*(x)^(a-1)*exp(-0.5*((lbar)^a+fi^a-(x)^a-fj^a)^2)*uk-c;

end

