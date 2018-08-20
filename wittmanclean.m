clear;
clc;
options = optimset('Display','off');

gamma = 1;       %Wittman specification parameter. 1 = Alessina, 0 = Standard
win = 0;            %Wittman win value
xihat = 0.3;
xjhat = 0.6;

t = 1000;

alpha = 2;
beta = 2;
weight = 1;

x1 = 0.01;

betamedian = @(x)betamedian(x,alpha,beta);
m = fsolve(betamedian,x1,options);

prob = zeros(1,t+1);
val = zeros(1,t+1);
pos = zeros(1,t+1);
brfi = zeros(1,t+1);
brfj = zeros(1,t+1);

for o = 1:t+1
    
    pos(o) = (o-1)/t;
    xj = (o-1)/t;
    
    for n = 1:t+1
       
       
        
        xi = (n-1)/t;
       
       
                      
       
        
        
        if xi<xj
            prob(n)=betacdf((xi+xj)/2,alpha,beta);
        elseif xi>xj
            prob(n)=(1-betacdf((xi+xj)/2,alpha,beta));
        else
            prob(n)=0.5;
        end
        
        val(n)=gamma*((prob(n)*(-(xi-xihat)^2+win)+(1-prob(n))*(-(xj-xihat)^2)))+(1-gamma)*prob(n);

        
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
       
        
        
        xj = (n-1)/t;
        
       
                      
       
        
        
        if xi<xj
            prob(n)=betacdf((xi+xj)/2,alpha,beta);
        elseif xi>xj
            prob(n)=(1-betacdf((xi+xj)/2,alpha,beta));
        else
            prob(n)=0.5;
        end
        
        
        val(n)=gamma*(prob(n)*(-(xi-xjhat)^2)+(1-prob(n))*(-(xj-xjhat)^2+win))+(1-gamma)*(1-prob(n));        
        
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

%filename = 'brf.xlsx';
%xlswrite(filename,transpose(pos),1,'A1')
%xlswrite(filename,transpose(brfi),1,'B1')
%xlswrite(filename,transpose(brfj),1,'C1')