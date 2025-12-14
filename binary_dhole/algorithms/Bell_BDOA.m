
function [Best_fitness,sFeat,Sf,Nf,cuve_f]  =Bell_BDOA(feat,label,N,max_Iter,alpha,HO)
%% Define Parameters
% Parameters
lb    = 0;
ub    = 1; 
thres = 0.5; 
fobj=@jFitnessFunction;
dim= size(feat,2);
cuve_f=zeros(1,N); 
X=initialization(N,dim,ub,lb); %Initialize population
Sol=zeros(N,dim);
global_Cov = zeros(1,max_Iter);
Best_fitness = inf;
best_position = zeros(1,dim);
fitness_f = zeros(1,N);

for i=1:N
   fitness_f(i) =  fobj(feat,label,(X(i,:) > thres),HO); %Calculate the fitness value of the function
   
   if fitness_f(i)<Best_fitness
       Best_fitness = fitness_f(i);
       localBest_position = X(i,:);
   end
end
global_position = localBest_position; 
global_fitness = Best_fitness;
cuve_f(1)=Best_fitness;
t=1; 
while(t<=max_Iter)
    C = 1-(t/max_Iter); %Eq.(7)
    PMN = round(rand*15+5); %Eq.(3)
    prey = (global_position+localBest_position)/2; %Eq.(5)
    prey_local = localBest_position;
        
    for i = 1:N
        
        if rand()<0.5
            if PMN<10
            %% Searching stage
              for j = 1:dim  
                Xnew(i,:) = X(i,:)+C*rand.*(prey(j)-X(i,:)); %Eq.(6)
              end
            else
            %% Encircling stage
                for j = 1:dim
                    z = round(rand*(N-1))+1;  %Eq.(9)
                    Xnew(i,j) = X(i,j)-X(z,j)+prey(j);  %Eq.(8)
                end
            end
        else
            %% Hunting stage
            %D_prey=global_position; %Eq.(10)
            Q = 3*rand*fitness_f(i)/fobj(feat,label,(prey_local > thres),HO); %Eq.(10)
            if Q>2   % The food is too big
                 W_prey = round(exp(-1/Q).*prey_local); %round is added   %Eq.(11)
                for j = 1:dim
                    Xnew(i,j) = X(i,j)+cos(2*pi*rand)*W_prey(j)*p_obj(PWN)-sin(2*pi*rand)*W_prey(j)*p_obj(PWN); %Eq.(12)
                end
            else
                Xnew(i,:) = (X(i,:)-prey_local)*p_obj(PWN)+p_obj(PWN).*rand(1,dim).*X(i,:); %Eq.(13)
            end
        end
    end
    %% boundary conditions
    for i=1:N
        for j =1:dim
            if length(ub)==1
                Xnew(i,j) = min(ub,Xnew(i,j));
                Xnew(i,j) = max(lb,Xnew(i,j));
            else
                Xnew(i,j) = min(ub(j),Xnew(i,j));
                Xnew(i,j) = max(lb(j),Xnew(i,j));
            end
        end
    end
    %% Binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xnew(1,:);%%local
    for i = 1:N
        for j= 1:dim
            Bell_shaped_transfer= abs(exp(-((Xnew(i,j))-1)^2)); %% Bell-shape, 1th function
            brandom=rand();
            if brandom > Bell_shaped_transfer
                Sol(i,j)=~Sol(i,j);
            else
                Sol(i,j)=Sol(i,j);
            end
        end
    end
    
    
    %% Updating new generation %%%%%%%%%%%%%
    localBest_position = Sol(1,:);%%local
    localBest_fitness = fobj(feat,label,(localBest_position > thres),HO);%%local
 
    for i =1:N
         %% Obtain the optimal solution for the updated population
        local_fitness = fobj(feat,label,(Sol(i,:)  > thres) ,HO);
        if local_fitness<localBest_fitness
                 localBest_fitness = local_fitness;
                 localBest_position = Sol(i,:);
        end
        %% Update the population to a new location
        if local_fitness<fitness_f(i)
             fitness_f(i) = local_fitness;
             X(i,:) = Sol(i,:);
             if fitness_f(i)<Best_fitness
                 Best_fitness=fitness_f(i);
                 global_position = X(i,:);
             end
        end
    end

    global_Cov(t) = localBest_fitness;
    cuve_f(t) = Best_fitness;
    t=t+1;
%     if mod(t,1)==0
%       disp("Bell_BDOA"+"iter"+num2str(t)+": "+Best_fitness); 
%    end
end
 best_fun = Best_fitness;
 % Selects features
Pos   = 1:dim;
Sf    = Pos((global_position > thres) == 1); 
sFeat = feat(:,Sf);
Nf    = length(Sf);
  fprintf('BDOA completed, ')
end
function y = p_obj(x)   %Eq.(4)
    PMN=x;
    C1 = 1;
    mu = 25;
    k = 0.5;
    D=rand;
    y = ((C1 / (1 + exp(-k * (PMN- mu))))^2)* D;
end



