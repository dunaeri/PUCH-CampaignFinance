clear;
clc;
options = optimset('Display','off');
format long;

t = 10;           %Number of iterations

gamma = 1;          %Wittman specification parameter. 1 = Wittman, 0 = Downs
win = 0;            %Wittman win value
xihat = 0.3;        %Candidate i's ideal position
xjhat = 0.7;        %Candidate j's ideal position
xk = 0.2;           %Donor k's ideal position 
xl = 0.8;           %Donor l's ideal position
kbar = 10;          %Donor k's contribution limit
lbar = 10;          %Donor l's contribution limit
eta = 0.5;          %Concavity parameter of donation function
c = 0.03;           %Marginal cost of donations
fi = 0;             %Public funding for candidate i
fj = 0;             %Public funding for candidate j

alpha = 2;          %Alpha parameter of beta distribution
beta = 2;           %Beta parameter of beta distribution
lambda = 0.5;       %Weight of donations on winning the election (NOTE: Inverse of lambda in paper. 1 = Donations not effective)

x1 = 0.01;          %Numerical solver parameters
x2 = [0.01,0.01];

prob = zeros(1,t+1);
val = zeros(1,t+1);
pos = zeros(1,t+1);
brfi = zeros(1,t+1);
brfj = zeros(1,t+1);

for o = 1:t+1
    
    pos(o) = (o-1)/t;
    xj = (o-1)/t;
    
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
    
    mx=1;
    for n = 1:t+1
        if val(n)>val(mx)
            mx=n;
        end
    end
    
    brfi(o)=(mx-1)/t;
    
    xi = (o-1)/t;
   
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
    
    mx=1;
    for n = 1:t+1
        if val(n)>val(mx)
            mx=n;
        end
    end
    
    brfj(o)=(mx-1)/t;
    
    brf = [pos(o),brfi(o),brfj(o)];
    disp(brf);
    
end

intersection = [0,0];

for o = 1:t+1
    if intersection ~= [0,0]
        break;
    end
    for n = 1:t+1
        if (pos(o)==brfj(n)) && (pos(n)==brfi(o))
            intersection = [brfi(o),brfj(n)];
            break;
        end
    end
end

disp(intersection);

filename = 'brf.xlsx';
xlswrite(filename,transpose(pos),1,'A1')
xlswrite(filename,transpose(brfi),1,'B1')
xlswrite(filename,transpose(brfj),1,'C1')