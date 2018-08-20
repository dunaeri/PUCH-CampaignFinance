clear;
clc;
options = optimset('Display','off');

gamma = 1;       %Wittman specification parameter. 1 = Wittman, 0 = Standard
win = 0;            %Wittman win value
xihat = 0.3;
xjhat = 0.7;

t = 1000;
u = 5;
xk = 0.2;
xl = 0.8;
kbar = 100;
lbar = 100;
eta = 0.5;

fi = 0;
fj = 0;

alpha = 2;
beta = 2;

linit = 1.25;
cinit = 0.06;

x1 = 0.01;
x2 = [0.01,0.01];

betamedian = @(x)betamedian(x,alpha,beta);
m = fsolve(betamedian,x1,options);

meas = [xk,xl,m];
low = min(meas);
high = max(meas);

prob = zeros(1,t+1);
val = zeros(1,t+1);
pos = zeros(1,t+1);
brfi = zeros(1,t+1);
brfj = zeros(1,t+1);
output = zeros(u^2,5);

lowrangei = 0.34;
lowrangej = 0.56;
highrangei = 0.44;
highrangej = 0.66;
radius = 0;

for q = 1:u
    c = cinit - q*0.01;
for v = 1:u

   lambda = linit - v*0.25;
    
prob = zeros(1,t+1);
val = zeros(1,t+1);
pos = zeros(1,t+1);
brfi = zeros(1,t+1);
brfj = zeros(1,t+1);



for o = 1:t+1
    
    val = zeros(1,t+1);
    pos(o) = (o-1)/t;
    xj = (o-1)/t;
    
    if (xj >= lowrangej - radius) && (xj <= highrangej + radius)
    for n = 1:t+1
       
        clearvars interiorkilj interiorki interiorki2 interiorlj interiorlj2 interiorkjli interiorkj interiorkj2 interiorli interiorli2 interiordi interiordj
        
        xi = (n-1)/t;
        
        uk = (xj-xi)*(xi+xj-2*xk);
        ul = (xj-xi)*(xi+xj-2*xl);
        ki = 0;
        kj = 0;
        li = 0;
        lj = 0;
                      
        if (uk>0) && (ul<0)
            interiorkilj = @(x)interiorkilj(x,eta,uk,ul,c,fi,fj);
            donation = fsolve(interiorkilj,x2,options);            
                
            if imag(sum(donation(:)))~=0
                ki=0;
                kj=0;
                li=0;
                lj=0;
            elseif (donation(1)>kbar && donation(2)>lbar)
                ki=kbar;
                lj=lbar;
            elseif donation(1)>kbar
                interiorlj = @(x)interiorlj(x,eta,ul,kbar,c,fi,fj);
                ki=kbar;
                lj=fsolve(interiorlj,x1,options);
                if lj>lbar
                    lj=lbar;
                end
            elseif donation(2)>lbar
                interiorki = @(x)interiorki(x,eta,uk,lbar,c,fi,fj);
                ki=fsolve(interiorki,x1,options);
                if ki>kbar
                    ki=kbar;
                end
                lj=lbar;
            elseif donation(1)<0
                interiorlj2 = @(x)interiorlj2(x,eta,ul,c,fi,fj);
                lj=fsolve(interiorlj2,x1,options);
            elseif donation(2)<0
                interiorki2 = @(x)interiorki2(x,eta,uk,c,fi,fj);
                ki=fsolve(interiorki2,x1,options);
            else
                ki=donation(1);
                lj=donation(2);
            end
        elseif (uk<0) && (ul>0)         
            interiorkjli = @(x)interiorkjli(x,eta,uk,ul,c,fi,fj);
            donation = fsolve(interiorkjli,x2,options);
            
            if imag(sum(donation(:)))~=0
                ki=0;
                kj=0;
                li=0;
                lj=0;
            elseif (donation(2)>kbar && donation(1)>lbar)
                kj=kbar;
                li=lbar;
            elseif donation(1)>lbar
                interiorkj = @(x)interiorkj(x,eta,uk,lbar,c,fi,fj);
                kj=fsolve(interiorkj,x1,options);
                li=lbar;
            elseif donation(2)>kbar
                interiorli = @(x)interiorli(x,eta,ul,kbar,c,fi,fj);
                kj=kbar;
                li=fsolve(interiorli,x1,options);
            elseif donation(1)<0
                interiorkj2 = @(x)interiorkj2(x,eta,uk,c,fi,fj);
                kj=fsolve(interiorkj2,x1,options);
            elseif donation(2)<0
                interiorli2 = @(x)interiorli2(x,eta,ul,c,fi,fj);
                li=fsolve(interiorli2,x1,options);
            else
                kj=donation(2);
                li=donation(1);
            end
        elseif (uk>=0) && (ul>=0)
            if (uk==0) && (ul==0)
                ki=0;
                li=0;
            elseif uk==0
                interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    li=0;
                elseif donation < lbar
                    li=donation;
                else
                    li = lbar;
                end
            elseif ul==0
                interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    ki=0;
                elseif donation < kbar
                    ki=donation;
                else
                    ki = kbar;
                end
            elseif uk>=ul
                interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= kbar
                    ki=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                    donation = fsolve(interiordi,x1,options);
                
                    if donation - kbar < 0
                        ki=kbar;
                    elseif donation - kbar <= lbar
                        ki=kbar;
                        li=donation-kbar;
                    else
                        ki=kbar;
                        li=lbar;
                    end
                end
            else
                interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= lbar
                    li=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                    donation = fsolve(interiordi,x1,options);
                
                    if donation - lbar < 0
                        li=lbar;
                    elseif donation - lbar <= kbar
                        li=lbar;
                        ki=donation-lbar;
                    else
                        li=lbar;
                        ki=kbar;
                    end
                end    
            end            
       elseif (uk<=0) && (ul<=0)
            if (uk==0) && (ul==0)
                kj=0;
                lj=0;
            elseif uk==0
                interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    lj=0;
                elseif donation < lbar
                    lj=donation;
                else
                    lj = lbar;
                end
            elseif ul==0
                interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    kj=0;
                elseif donation < kbar
                    kj=donation;
                else
                    kj = kbar;
                end
            elseif uk<=ul
                interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= kbar
                    kj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                    donation = fsolve(interiordj,x1,options);
                
                    if donation - kbar < 0
                        kj=kbar;
                    elseif donation - kbar <= lbar
                        kj=kbar;
                        lj=donation-kbar;
                    else
                        kj=kbar;
                        lj=lbar;
                    end
                end
            else
                interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= lbar
                    lj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                    donation = fsolve(interiordj,x1,options);
                
                    if donation - lbar < 0
                        lj=lbar;
                    elseif donation - lbar <= kbar
                        lj=lbar;
                        kj=donation-lbar;
                    else
                        lj=lbar;
                        kj=kbar;
                    end
                end    
            end            
        end
        
       if xi<xj
            prob(n)=lambda*betacdf((xi+xj)/2,alpha,beta)+(1-lambda)*normcdf((ki+li)^eta+fi^eta-(kj+lj)^eta-fj^eta);
        elseif xi>xj
            prob(n)=lambda*(1-betacdf((xi+xj)/2,alpha,beta))+(1-lambda)*normcdf((ki+li)^eta+fi^eta-(kj+lj)^eta-fj^eta);
        else
            prob(n)=0.5;
        end
        
        val(n)=(prob(n)*(-(xi-xihat)^2+win)+(1-prob(n))*(-(xj-xihat)^2))*gamma + prob(n)*(1-gamma);

        
    end
    end
    
    mx=1;
    for n = 1:t+1
        if val(n)>val(mx)
            mx=n;
        end
    end
    
    brfi(o)=(mx-1)/t;
    
    val = zeros(1,t+1);
    xi = (o-1)/t;
   
    if (xi >= lowrangei - radius) && (xi <= highrangei + radius)
    for n = 1:t+1
       
        clearvars interiorkilj interiorki interiorki2 interiorlj interiorlj2 interiorkjli interiorkj interiorkj2 interiorli interiorli2 interiordi interiordj
        
        xj = (n-1)/t;
        uk = (xj-xi)*(xi+xj-2*xk);
        ul = (xj-xi)*(xi+xj-2*xl);
        ki = 0;
        kj = 0;
        li = 0;
        lj = 0;
                      
        if (uk>0) && (ul<0)
            interiorkilj = @(x)interiorkilj(x,eta,uk,ul,c,fi,fj);
            donation = fsolve(interiorkilj,x2,options);            
                
            if imag(sum(donation(:)))~=0
                ki=0;
                kj=0;
                li=0;
                lj=0;
            elseif (donation(1)>kbar && donation(2)>lbar)
                ki=kbar;
                lj=lbar;
            elseif donation(1)>kbar
                interiorlj = @(x)interiorlj(x,eta,ul,kbar,c,fi,fj);
                ki=kbar;
                lj=fsolve(interiorlj,x1,options);
                if lj>lbar
                    lj=lbar;
                end
            elseif donation(2)>lbar
                interiorki = @(x)interiorki(x,eta,uk,lbar,c,fi,fj);
                ki=fsolve(interiorki,x1,options);
                if ki>kbar
                    ki=kbar;
                end
                lj=lbar;
            elseif donation(1)<0
                interiorlj2 = @(x)interiorlj2(x,eta,ul,c,fi,fj);
                lj=fsolve(interiorlj2,x1,options);
            elseif donation(2)<0
                interiorki2 = @(x)interiorki2(x,eta,uk,c,fi,fj);
                ki=fsolve(interiorki2,x1,options);
            else
                ki=donation(1);
                lj=donation(2);
            end
        elseif (uk<0) && (ul>0)         
            interiorkjli = @(x)interiorkjli(x,eta,uk,ul,c,fi,fj);
            donation = fsolve(interiorkjli,x2,options);
            
            if imag(sum(donation(:)))~=0
                ki=0;
                kj=0;
                li=0;
                lj=0;
            elseif (donation(2)>kbar && donation(1)>lbar)
                kj=kbar;
                li=lbar;
            elseif donation(1)>lbar
                interiorkj = @(x)interiorkj(x,eta,uk,lbar,c,fi,fj);
                kj=fsolve(interiorkj,x1,options);
                li=lbar;
            elseif donation(2)>kbar
                interiorli = @(x)interiorli(x,eta,ul,kbar,c,fi,fj);
                kj=kbar;
                li=fsolve(interiorli,x1,options);
            elseif donation(1)<0
                interiorkj2 = @(x)interiorkj2(x,eta,uk,c,fi,fj);
                kj=fsolve(interiorkj2,x1,options);
            elseif donation(2)<0
                interiorli2 = @(x)interiorli2(x,eta,ul,c,fi,fj);
                li=fsolve(interiorli2,x1,options);
            else
                kj=donation(2);
                li=donation(1);
            end
        elseif (uk>=0) && (ul>=0)
            if (uk==0) && (ul==0)
                ki=0;
                li=0;
            elseif uk==0
                interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    li=0;
                elseif donation < lbar
                    li=donation;
                else
                    li = lbar;
                end
            elseif ul==0
                interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    ki=0;
                elseif donation < kbar
                    ki=donation;
                else
                    ki = kbar;
                end
            elseif uk>=ul
                interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= kbar
                    ki=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                    donation = fsolve(interiordi,x1,options);
                
                    if donation - kbar < 0
                        ki=kbar;
                    elseif donation - kbar <= lbar
                        ki=kbar;
                        li=donation-kbar;
                    else
                        ki=kbar;
                        li=lbar;
                    end
                end
            else
                interiordi = @(x)interiordi(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= lbar
                    li=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,uk,c,fi,fj);
                    donation = fsolve(interiordi,x1,options);
                
                    if donation - lbar < 0
                        li=lbar;
                    elseif donation - lbar <= kbar
                        li=lbar;
                        ki=donation-lbar;
                    else
                        li=lbar;
                        ki=kbar;
                    end
                end    
            end            
       elseif (uk<=0) && (ul<=0)
            if (uk==0) && (ul==0)
                kj=0;
                lj=0;
            elseif uk==0
                interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    lj=0;
                elseif donation < lbar
                    lj=donation;
                else
                    lj = lbar;
                end
            elseif ul==0
                interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    kj=0;
                elseif donation < kbar
                    kj=donation;
                else
                    kj = kbar;
                end
            elseif uk<=ul
                interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= kbar
                    kj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                    donation = fsolve(interiordj,x1,options);
                
                    if donation - kbar < 0
                        kj=kbar;
                    elseif donation - kbar <= lbar
                        kj=kbar;
                        lj=donation-kbar;
                    else
                        kj=kbar;
                        lj=lbar;
                    end
                end
            else
                interiordj = @(x)interiordj(x,eta,ul,c,fi,fj);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= lbar
                    lj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,uk,c,fi,fj);
                    donation = fsolve(interiordj,x1,options);
                
                    if donation - lbar < 0
                        lj=lbar;
                    elseif donation - lbar <= kbar
                        lj=lbar;
                        kj=donation-lbar;
                    else
                        lj=lbar;
                        kj=kbar;
                    end
                end    
            end            
        end
        
        if xi<xj
            prob(n)=lambda*betacdf((xi+xj)/2,alpha,beta)+(1-lambda)*normcdf((ki+li)^eta+fi^eta-(kj+lj)^eta-fj^eta);
        elseif xi>xj
            prob(n)=lambda*(1-betacdf((xi+xj)/2,alpha,beta))+(1-lambda)*normcdf((ki+li)^eta+fi^eta-(kj+lj)^eta-fj^eta);
        else
            prob(n)=0.5;
        end
        
        val(n)=(prob(n)*(-(xi-xjhat)^2)+(1-prob(n))*(-(xj-xjhat)^2)+win)*gamma + (1-prob(n))*(1-gamma);        
        
    end
    end
    
    mx=1;
    for n = 1:t+1
        if val(n)>val(mx)
            mx=n;
        end
    end
    
    brfj(o)=(mx-1)/t;
    
    brf = [pos(o),brfi(o),brfj(o)];
    
end

intersection = [0,0];

for o = 1:t+1
    if intersection ~= [0,0]
        break;
    end
    for n = 1:t+1
        if (pos(o)==brfj(n)) && (pos(n)==brfi(o)) && (pos(o)~=0) && (pos(n)~=0)
            intersection = [brfi(o),brfj(n)];
            break;
        end
    end
end

outputnum = 5*(q-1)+v;

output(outputnum,1) = c;
output(outputnum,2) = 1 - lambda;
output(outputnum,3) = intersection(1);
output(outputnum,4) = intersection(2);
output(outputnum,5) = intersection(2)-intersection(1);
disp(output(outputnum,:));

end
end

filename = 'divergence.xlsx';
xlswrite(filename,output,1,'A1')

