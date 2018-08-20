clear;
clc;
options = optimset('Display','off');

xiinit = 0.686;    
xjinit = 0.686; 

gamma = 0;       %Wittman specification parameter. 1 = Wittman, 0 = Downs
win = 0;            %Wittman win value
xihat = 0.3;
xjhat = 0.7;

t = 1000;
xk = 0.2;
xl = 0.8;
kbar = 1000000;
lbar = 1000000;
eta = 0.5;
c = 0.06;           %Marginal cost of donations

alpha = 4;
beta = 2;
lambda = 0.5;

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
brfi = zeros(6,t+1);
brfj = zeros(6,t+1);
mx = zeros(1,2);   

if xjinit ~= 0
    
    xj = xjinit;
    
    for n = 1:t+1
       
        clearvars interiorkilj interiorki interiorki2 interiorlj interiorlj2 interiorkjli interiorkj interiorkj2 interiorli interiorli2 interiordi interiordj
        
        pos(n) = (n-1)/t;
        
        xi = (n-1)/t;
        uk = (xj-xi)*(xi+xj-2*xk);
        ul = (xj-xi)*(xi+xj-2*xl);
        ki = 0;
        kj = 0;
        li = 0;
        lj = 0;
                      
        if (uk>0) && (ul<0)
            interiorkilj = @(x)interiorkilj(x,eta,uk,ul,c);
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
                interiorlj = @(x)interiorlj(x,eta,ul,kbar,c);
                ki=kbar;
                lj=fsolve(interiorlj,x1,options);
                if lj>lbar
                    lj=lbar;
                end
            elseif donation(2)>lbar
                interiorki = @(x)interiorki(x,eta,uk,lbar,c);
                ki=fsolve(interiorki,x1,options);
                if ki>kbar
                    ki=kbar;
                end
                lj=lbar;
            elseif donation(1)<0
                interiorlj2 = @(x)interiorlj2(x,eta,ul,c);
                lj=fsolve(interiorlj2,x1,options);
            elseif donation(2)<0
                interiorki2 = @(x)interiorki2(x,eta,uk,c);
                ki=fsolve(interiorki2,x1,options);
            else
                ki=donation(1);
                lj=donation(2);
            end
        elseif (uk<0) && (ul>0)         
            interiorkjli = @(x)interiorkjli(x,eta,uk,ul,c);
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
                interiorkj = @(x)interiorkj(x,eta,uk,lbar,c);
                kj=fsolve(interiorkj,x1,options);
                li=lbar;
            elseif donation(2)>kbar
                interiorli = @(x)interiorli(x,eta,ul,kbar,c);
                kj=kbar;
                li=fsolve(interiorli,x1,options);
            elseif donation(1)<0
                interiorkj2 = @(x)interiorkj2(x,eta,uk,c);
                kj=fsolve(interiorkj2,x1,options);
            elseif donation(2)<0
                interiorli2 = @(x)interiorli2(x,eta,ul,c);
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
                interiordi = @(x)interiordi(x,eta,ul,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    li=0;
                elseif donation < lbar
                    li=donation;
                else
                    li = lbar;
                end
            elseif ul==0
                interiordi = @(x)interiordi(x,eta,uk,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    ki=0;
                elseif donation < kbar
                    ki=donation;
                else
                    ki = kbar;
                end
            elseif uk>=ul
                interiordi = @(x)interiordi(x,eta,uk,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= kbar
                    ki=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,ul,c);
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
                interiordi = @(x)interiordi(x,eta,ul,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= lbar
                    li=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,uk,c);
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
                interiordj = @(x)interiordj(x,eta,ul,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    lj=0;
                elseif donation < lbar
                    lj=donation;
                else
                    lj = lbar;
                end
            elseif ul==0
                interiordj = @(x)interiordj(x,eta,uk,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    kj=0;
                elseif donation < kbar
                    kj=donation;
                else
                    kj = kbar;
                end
            elseif uk<=ul
                interiordj = @(x)interiordj(x,eta,uk,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= kbar
                    kj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,ul,c);
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
                interiordj = @(x)interiordj(x,eta,ul,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= lbar
                    lj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,uk,c);
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
            prob(n)=lambda*betacdf((xi+xj)/2,alpha,beta)+(1-lambda)*normcdf((ki+li)^eta-(kj+lj)^eta);
        elseif xi>xj
            prob(n)=lambda*(1-betacdf((xi+xj)/2,alpha,beta))+(1-lambda)*normcdf((ki+li)^eta-(kj+lj)^eta);
        else
            prob(n)=0.5;
        end
        
        val(n)=(prob(n)*(-(xi-xihat)^2+win)+(1-prob(n))*(-(xj-xihat)^2))*gamma + prob(n)*(1-gamma);

        brfi(1,n)=val(n);
        brfi(2,n)=ki;
        brfi(3,n)=li;
        brfi(4,n)=kj;
        brfi(5,n)=lj;
        brfi(6,n)=prob(n);
        
    end
   
    mx(1,1)=1;
    for n = 1:t+1
        if val(n)>val(mx(1,1))
            mx(1,1)=n;
        end
    end
    mx(1,1)=(mx(1,1)-1)/t;
    
end
    
 if xiinit ~= 0   
    xi = xiinit;
   
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
            interiorkilj = @(x)interiorkilj(x,eta,uk,ul,c);
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
                interiorlj = @(x)interiorlj(x,eta,ul,kbar,c);
                ki=kbar;
                lj=fsolve(interiorlj,x1,options);
                if lj>lbar
                    lj=lbar;
                end
            elseif donation(2)>lbar
                interiorki = @(x)interiorki(x,eta,uk,lbar,c);
                ki=fsolve(interiorki,x1,options);
                if ki>kbar
                    ki=kbar;
                end
                lj=lbar;
            elseif donation(1)<0
                interiorlj2 = @(x)interiorlj2(x,eta,ul,c);
                lj=fsolve(interiorlj2,x1,options);
            elseif donation(2)<0
                interiorki2 = @(x)interiorki2(x,eta,uk,c);
                ki=fsolve(interiorki2,x1,options);
            else
                ki=donation(1);
                lj=donation(2);
            end
        elseif (uk<0) && (ul>0)         
            interiorkjli = @(x)interiorkjli(x,eta,uk,ul,c);
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
                interiorkj = @(x)interiorkj(x,eta,uk,lbar,c);
                kj=fsolve(interiorkj,x1,options);
                li=lbar;
            elseif donation(2)>kbar
                interiorli = @(x)interiorli(x,eta,ul,kbar,c);
                kj=kbar;
                li=fsolve(interiorli,x1,options);
            elseif donation(1)<0
                interiorkj2 = @(x)interiorkj2(x,eta,uk,c);
                kj=fsolve(interiorkj2,x1,options);
            elseif donation(2)<0
                interiorli2 = @(x)interiorli2(x,eta,ul,c);
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
                interiordi = @(x)interiordi(x,eta,ul,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    li=0;
                elseif donation < lbar
                    li=donation;
                else
                    li = lbar;
                end
            elseif ul==0
                interiordi = @(x)interiordi(x,eta,uk,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation < 0
                    ki=0;
                elseif donation < kbar
                    ki=donation;
                else
                    ki = kbar;
                end
            elseif uk>=ul
                interiordi = @(x)interiordi(x,eta,uk,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= kbar
                    ki=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,ul,c);
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
                interiordi = @(x)interiordi(x,eta,ul,c);
                donation = fsolve(interiordi,x1,options);
            
                if donation <= lbar
                    li=donation;
                else
                    clearvars interiordi
                    interiordi = @(x)interiordi(x,eta,uk,c);
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
                interiordj = @(x)interiordj(x,eta,ul,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    lj=0;
                elseif donation < lbar
                    lj=donation;
                else
                    lj = lbar;
                end
            elseif ul==0
                interiordj = @(x)interiordj(x,eta,uk,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation < 0
                    kj=0;
                elseif donation < kbar
                    kj=donation;
                else
                    kj = kbar;
                end
            elseif uk<=ul
                interiordj = @(x)interiordj(x,eta,uk,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= kbar
                    kj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,ul,c);
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
                interiordj = @(x)interiordj(x,eta,ul,c);
                donation = fsolve(interiordj,x1,options);
            
                if donation <= lbar
                    lj=donation;
                else
                    clearvars interiordj
                    interiordj = @(x)interiordj(x,eta,uk,c);
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
            prob(n)=lambda*betacdf((xi+xj)/2,alpha,beta)+(1-lambda)*normcdf((ki+li)^eta-(kj+lj)^eta);
        elseif xi>xj
            prob(n)=lambda*(1-betacdf((xi+xj)/2,alpha,beta))+(1-lambda)*normcdf((ki+li)^eta-(kj+lj)^eta);
        else
            prob(n)=0.5;
        end
        
        val(n)=(prob(n)*(-(xi-xjhat)^2)+(1-prob(n))*(-(xj-xjhat)^2)+win)*gamma + (1-prob(n))*(1-gamma);        
        
        brfj(1,n)=val(n);
        brfj(2,n)=ki;
        brfj(3,n)=li;
        brfj(4,n)=kj;
        brfj(5,n)=lj;
        brfj(6,n)=prob(n);

        
    end
    
 mx(1,2)=1;
    for n = 1:t+1
        if val(n)>val(mx(1,2))
            mx(1,2)=n;
        end
    end
 mx(1,2)=(mx(1,2)-1)/t;
    
 end
    
 disp(mx);   
    
    

filename = 'brfzoom.xlsx';
xlswrite(filename,transpose(pos),1,'A2')
xlswrite(filename,transpose(brfi),1,'B2')
xlswrite(filename,transpose(pos),2,'A2')
xlswrite(filename,transpose(brfj),2,'B2')