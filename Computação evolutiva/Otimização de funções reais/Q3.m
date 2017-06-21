close all;
clc;
clear;

filename = 'exe1.mat';
n = 100; % número de indivíduos
cr = 0.8; % probabilidade de crossover 
m = 0.2; % probabilidade de mutação
maxIt = 1000; % número máximo de iterações
q=2; %número de indiníduos no torneio
sigma = 0.5;
alpha = 1;
gamaMate = 0.5;

mFit =0; % Melhor valor de fitness
mIndv = [0,0];

x =  -1 + 3* rand(1,n);
y =  -1 + 3* rand(1,n);

p = [x' y'];



if(nnz(p(:,2)<-1)>0)
    disp('error inicialização');
    break;
end

for i=1:maxIt
     i
    %crossOver ----------------------
    pc = rand(n,1);
    crossIdx = find(pc>(1-cr));
    while(length(crossIdx)>0)
        if(length(crossIdx)==1)
            crossIdx(end+1) = randi(n,1); 
        end
        a = randperm(length(crossIdx),2);
        
        %DEB
        a3 = a(2);
        while (true)
            if (sqrt( (p(a(1),1)-p(a3,1))^2 + (p(a(1),2)-p(a3,2))^2) >gamaMate)
                a3 = randi(n,1);
            else
                break;
            end
        end
        % Recombinação intermediária local
        p(end+1,:) = [ (p(a(1),1)+p(a(2),1))/2 (p(a(1),2)+p(a(2),2))/2];
        
        
        
        %p(end+1) = [ (p(a(1),1)*0.6) + (p(a(2),1)*0.4) p(a(1),2)*0.6) + (p(a(2),2)*0.4)];
        %p(end+1) = [ (p(a(1),1)*0.4) + (p(a(2),1)*0.6) p(a(1),2)*0.4) + (p(a(2),2)*0.6)];
        
        crossIdx([a(1),a(2)])=[];
    end
    %-----------------------------------
    
    if(nnz(p(:,2)<-1)>0)
    disp('error cross');
    break;
    end
    
     %Mutação------------------------------------ 
     pm = rand(n,1);
     mutacaoIdx = find(pm>(1-m));
     
     while(length(mutacaoIdx)>0)
        mf = normrnd(0,1,1,2);
        a = randperm(length(mutacaoIdx),1);
        p(a,1) = p(a,1) + mf(1);
        p(a,2) = p(a,2) + mf(2);
        
        if(p(a,1)>2)
            p(a,1)=2;
        end
        if(p(a,2)>2)
            p(a,2)=2;
        end
        if(p(a,1)<-1)
            p(a,1)=-1;
        end
        if(p(a,2)<-1)
            p(a,2)=-1;
        end
        
        mutacaoIdx(a)=[];
     end
     %------------------------------------------
     if(nnz(p(:,2)<-1)>0)
        disp('error mutação');
        break;
     end
     
     %Fitness da população --------------------------------
     fp = zeros(size(p,1),1); % fp = fitness da população
     for j=1:size(p,1);
        fp(j)= p(j,1)*sin(4*pi*p(j,1)) - p(j,2)*sin(4*pi*p(j,2) + pi) +1;
     end
    %-----------------------------------------------------
    
    %Fitness sharing ----------------------------------
    
    nincho = zeros(size(p,1),1);
    for j=1:size(p,1);
        V = [ sqrt((p(:,1) - p(j,1)).^2  + (p(:,2) - p(j,2)).^2)];
        V(j) = sigma;
        D = V(find(V<sigma));
        nincho(j) = sum(1-(D/0.5))^(alpha);
    end
    %---------------------------------------------
    
    %rebalanceando fitness-----------------------
    fshp = fp ./ nincho;
    %------------------------------------------
    
    % Melhor individuo ---------------------
    if(max(fshp) > mFit)
        [mv, maxIdx]=max(fshp);
        mFit = mv;
        mIndv = p(maxIdx,:);
    end
    % ---------------------------------------
     
    %Seleção por torneio -------------------------
    [np,fnp]=SelecaoTorneioReal( q,p,fshp,n );
    % ------------------------------------------ 
    
    % Salvaciosnismo -----------------------------
    if (mFit>max(fnp))
        np(end+1,:)=mIndv;
        fnp(end+1)=mv;
    end
    % --------------------------------------------
    
    % médio, máximo e mínimo fitness ---------------
    meanP(i) = mean(fnp);
    maxP(i) = max(fnp);
    minP(i) = min(fnp);
    if(i>1)
     if(maxP(i)<maxP(i-1));
        disp('bugg');
        break;
     end
    end
    %----------------------------------------------    
    p = np;
end

%plotando gráficos------------------------
figure (1);
hold on;
plot(maxP, 'r-')
plot(meanP, 'g-')


%plotando superfície --------------------------------------------
X=-1:0.05:2;
Y=-1:0.05:2;

l = length(X);
MX = ones(l);
MY = ones(l);
MZ = ones(l);

for k=1:l
    MX(k,:) = X(k);
    MY(k,:) = Y(k);
end

MY = MY';

for k=1:l
    for t=1:l
        MZ(k,t) = MX(k,t)*sin(4*pi*MX(k,t)) - MY(k,t)*sin(4*pi*MY(k,t)+pi)+1;
    end
end

figure(2);
surf(MX,MY,MZ);
%---------------------------------------------------

% plotando indivíduos finais ----------------------
for j=1:size(p,1)
    hold on;
    plot3(p(j,1),p(j,2),p(j,1)*sin(4*pi*p(j,1)) - p(j,2)*sin(4*pi*p(j,2)+pi)+1, 'r*');
end    

% ----------------------------------------------



