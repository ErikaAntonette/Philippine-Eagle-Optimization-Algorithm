function [fbest_pheagle, xbest_pheagle, evalnum_pheagle] = pheaglealgorithm(D, f, Space_x_max, Space_x_min,IES,IFS,MFE)

%close all;

tic

%% Starting Parameters
if nargin < 7
    InitEagleSize = 20*D^2;  
    FoodSize = 10*D^2; 
%    MaxEvals = 1000*D;
     MaxEvals = 10000*D;
else
    InitEagleSize = IES;  
    FoodSize = IFS; 
    MaxEvals = MFE;
end

ClusterSize = max(0.02*min(Space_x_max - Space_x_min),1);
PS1 = InitEagleSize;
MinPopSize = 5;

%% Generation Number, Number of Function Evaluations
G = 0; 
CurrentEvals = 0;

%% Probability of Each Operator
Num_Operators = 3;
Prob_Operators = 1./Num_Operators.*ones(1,Num_Operators);

%% Archive Parameters
arch_rate = 2.6;
archive.NP = round(arch_rate * PS1); % the maximum size of the archive
archive.pop = zeros(0, D); % the solutions stored in the archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions

%% Plot Fitness Function
% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
% x = linspace(Space_x_min(1), Space_x_max(1), 200);
% y = linspace(Space_x_min(2), Space_x_max(2), 200);
% for i = 1 : size(x,2)
%     for j = 1: size(y,2)
%         z(j,i) = f([x(i),y(j)]');
%     end
% end
% figure(1)
% surfc(x,y,z)
% hold on;
% scatter3(f_optimalminimizer(1),f_optimalminimizer(2),f_optimalvalue,70,'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])

%% ===== START OF PHILIPPINE EAGLE ALGORITHM ===== %%

%% Initial Eagles
Eagle = repmat(Space_x_min,InitEagleSize,1)+repmat((Space_x_max-Space_x_min),InitEagleSize,1).*lhsdesign(InitEagleSize, D);

%% Calculate Fitness of Initial Eagles
FitEagle = zeros(InitEagleSize, 1);
for EagleNum = 1 : InitEagleSize
        FitEagle(EagleNum) = f(Eagle(EagleNum,:));
end
CurrentEvals = CurrentEvals + InitEagleSize;

%% Store Best Initial Eagle
[FitEagle, Index] = sort(FitEagle);
Eagle = Eagle(Index,:);
Eagle_Star = Eagle(1,:);
bestold = FitEagle(1);

%% Plot Initial Eagles, Best Eagle
% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
% figure(2);
% contour(x,y,z,35,'LevelStep',1)
% hold on;
% scatter3(f_optimalminimizer(1),f_optimalminimizer(2),f_optimalvalue,70,'filled','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% hold on; E = plot(InitEagle(:,1),InitEagle(:,2),'.','markersize',15); 
% pause(0.1);

%% Food Search of Best Eagle
Cluster_x_min = Eagle_Star-repmat(ClusterSize,1,D);
Cluster_x_max = Eagle_Star+repmat(ClusterSize,1,D);

options = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunctionEvaluations',FoodSize,'Display','off');
[Food_Star,FitFood] = fmincon(f,Eagle_Star,[],[],[],[],max(Cluster_x_min,Space_x_min), min(Cluster_x_max,Space_x_max),[],options);

CurrentEvals = CurrentEvals + FoodSize;

%% Choose Minimum Between Best Eagle and Best Food
%[BestValue, ~] = sort([FitEagle(1);FitFood]);
% x = [Eagle_Star;Food_Star];
% Star = x(BestIndex(1),:);

%% Plot Best Eagle/Food
% hold on; ES = plot(Eagle_Star(1), Eagle_Star(2),'h','markersize',13,'LineWidth',1.5);
% hold on; S = plot(Star(1), Star(2),'*','markersize',13,'LineWidth',1.5); pause(0.1);

%% Record Best Eagle and Initial Best Food
bestcost = FitFood;
Result = [];
Result = [Result; G Eagle_Star bestold Food_Star FitFood bestcost CurrentEvals];

%% For Adaptation of Scaling Factor F
hist_pos = 1;
memory_size = 20*D;
archive_f = ones(1,memory_size).*0.2;

%% Eagle Matrix Placeholder
% Eagle = Eagle; %to store initial eagles

%% ========================= START OF MAIN LOOP ========================= %%
%for GenNum = 1:MaxG
while CurrentEvals < MaxEvals

%     if abs(bestcost - f_optimalvalue) < 1e-08
%         break
%     end
    
G = G + 1;

%% Linear Reduction of PS1 
UpdPopSize = floor((((MinPopSize - InitEagleSize) / MaxEvals) * CurrentEvals) + InitEagleSize);
if PS1 > UpdPopSize
    reduction_ind_num = PS1 - UpdPopSize;
    if PS1 - reduction_ind_num <  MinPopSize
        reduction_ind_num = PS1 - MinPopSize;
    end
    %% remove the worst ind.
    for r = 1 : reduction_ind_num
        vv=PS1;
        Eagle(vv,:)=[];
        FitEagle(vv)=[];
        PS1 = PS1 - 1;
    end
    archive.NP = round(arch_rate * PS1);
    if size(archive.pop, 1) > archive.NP
        rndpos = randperm(size(archive.pop, 1));
        rndpos = rndpos(1 : archive.NP);
        archive.pop = archive.pop(rndpos, :);
    end
end

%% Calculate Scaling Factor F
mem_rand_index = ceil(memory_size * rand(PS1, 1));
mu_sf = archive_f(mem_rand_index);

F = mu_sf + 0.1 * tan(pi * (rand(1,PS1) - 0.5)); %Cauchy distribution, mean mu_sf, variance 0.1
pos = find(F <= 0);

while ~ isempty(pos)
    F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(1,length(pos)) - 0.5));
    pos = find(F <= 0);
end

F = min(F, 1);
F = F';

%% ================================= APPLICATION OF OPERATORS =================================
bb = rand(PS1, 1);
probiter = Prob_Operators(1,:);
l2 = sum(Prob_Operators(1:2));
op_1 = bb <=  probiter(1)*ones(PS1, 1);
op_2 = bb > probiter(1)*ones(PS1, 1) &  bb <= (l2*ones(PS1, 1)) ;
op_3 = bb > l2*ones(PS1, 1) &  bb <= (ones(PS1, 1)) ;

%% Compute Relative Distance of Each Eagle
RelDistance = squareform(pdist(Eagle));
RelDistance(RelDistance == 0) = NaN;
    
%% Assign Closest Eagle to Each Eagle
[~,Ind] = min(RelDistance);
EaglePairs = [1:PS1;Ind]';

%% Archive Initial Eagles
popAll = [Eagle;archive.pop];

%% Generate Random Indexing Numbers of Distinct Permutation
r0 = 1 : PS1;
r4 = EaglePairs(:,2)';

[r1, r2, r3] = gnR1R2_v2(PS1, size(popAll, 1), r0, r4);
    
%% ========== Movement Operator ========== %%
% DE current-to-phi_best-with-archive Operator for Movement of Eagles
EagleOffspring = zeros(PS1,D);
Eagle_Close = Eagle(EaglePairs(:,2),:);
beta = 1*exp(-1*(sqrt(sum((Eagle_Close - Eagle).^2,2))).^2);
Eagle_phi = repmat(Eagle_Star,PS1,1);
% X = 1-2*F;
% X = min(X,1);
% X = max(X,0);
        
EagleOffspring(op_1==1,:) = Eagle(op_1==1,:)... 
                            + F(op_1==1, ones(1, D)).*(Eagle_phi(op_1==1,:) - Eagle(op_1==1,:) + Eagle(r1(op_1==1),:) - popAll(r2(op_1==1),:))...
                            + F(op_1==1, ones(1, D)).*beta(op_1==1, ones(1, D)).*(Eagle_Close(op_1==1,:) - Eagle(op_1==1,:));

%% ========== IMODE Mutation Operator ========== %%
%EagleOffspring(op_2==1,:) = F(op_2==1, ones(1, D)).*Eagle(r1(op_2==1),:) + F(op_2==1, ones(1, D)).*(Eagle_phi(op_2==1,:) - Eagle(r3(op_2==1), :));
levyflightterm = rand(PS1,D).*repmat(Levy(D),PS1,1);

EagleOffspring(op_2==1,:) = F(op_2==1, ones(1, D)).*Eagle(r1(op_2==1),:) + F(op_2==1, ones(1, D)).*(Eagle_phi(op_2==1,:) - Eagle(r3(op_2==1), :)) + levyflightterm(op_2==1,:);
%.*(Eagle_phi(op_2==1,:) - Eagle(op_2==1,:));                    
      
%% ========== Harris Hawks Mutation Operator ========== %%
% rho = repmat(CurrentEvals/MaxEvals,PS1,1);
% Y = 1-rho;
Eagle_mean = repmat(mean(Eagle),PS1,1);
Z = repmat(Space_x_min,PS1,1)+repmat((Space_x_max-Space_x_min),PS1,1).*lhsdesign(PS1, D);

EagleOffspring(op_3==1,:) = F(op_3==1, ones(1, D)).*Z(op_3==1,:) + F(op_3==1, ones(1, D)).*(Eagle_phi(op_3==1,:) - Eagle_mean(op_3==1,:));

%% Handle Boundaries
EagleOffspring = han_boun(EagleOffspring, Space_x_max, Space_x_min, Eagle, PS1, 1);

%% Calculate Fitness of Eagle Offspring
NewFitEagle = zeros(PS1, 1);
for EagleNum = 1 : PS1
        NewFitEagle(EagleNum) = f(EagleOffspring(EagleNum,:));
end
CurrentEvals = CurrentEvals + PS1;

%% Calculate Improvement for Scaling Factor F
% bestold = FitEagle(1);
diff = abs(FitEagle - NewFitEagle);
I = (NewFitEagle < FitEagle);
goodF = F(I == 1);

%% Update Archive
archive = updateArchive(archive, Eagle(I == 1, :), FitEagle(I == 1));

%% Update Probability of Each Operator
diff2 = max(0,(FitEagle - NewFitEagle))./abs(FitEagle);
count_S(1) = max(0,mean(diff2(op_1==1)));
count_S(2) = max(0,mean(diff2(op_2==1)));
count_S(3) = max(0,mean(diff2(op_3==1)));

if count_S~=0
    Prob_Operators = max(0.1,min(0.9,count_S./(sum(count_S))));
else
    Prob_Operators = 1/3 * ones(1,3);
end

%% Update Eagles and Fitness of Eagles
FitEagle(I == 1) = NewFitEagle(I == 1);
Eagle(I == 1, :) = EagleOffspring(I == 1, :);
% NewEagle_Star = Eagle(1,:);

% %% Plot Eagle Offspring
% delete(E); 
% E = plot(Eagle(:,1),Eagle(:,2),'.','markersize',15);
% pause(0.1);

%% Update Memory of Scaling Factor F
if size(goodF,1)==1
    goodF = goodF';
end

num_success_params = numel(goodF);
if num_success_params > 0
    weightsDE = diff(I == 1)./ sum(diff(I == 1));    
    archive_f(hist_pos) = (weightsDE' * (goodF .^ 2))./ (weightsDE' * goodF);
    hist_pos = hist_pos + 1; 
    
    if hist_pos > memory_size  
        hist_pos = 1; 
    end
else
    archive_f(hist_pos)=0.5;
end

%% Sort New Eagle, New Fitness
[FitEagle, ind]=sort(FitEagle);
Eagle=Eagle(ind,:);

if FitEagle(1) < bestold
    NewEagle_Star = Eagle(1,:);
    Fitness = FitEagle(1);
else
    NewEagle_Star = Eagle_Star;
    Fitness = bestold;
end

%% Food Search of Best Eagle
if Fitness == bestold
    Cluster_x_min = Food_Star-repmat(ClusterSize,1,D);
    Cluster_x_max = Food_Star+repmat(ClusterSize,1,D);
    
    [New_Food_Star,NewFitFood] = fmincon(f,Food_Star,[],[],[],[],max(Cluster_x_min,Space_x_min), min(Cluster_x_max,Space_x_max),[],options);
else   
    Cluster_x_min = NewEagle_Star-repmat(ClusterSize,1,D);
    Cluster_x_max = NewEagle_Star+repmat(ClusterSize,1,D);

    [New_Food_Star,NewFitFood] = fmincon(f,NewEagle_Star,[],[],[],[],max(Cluster_x_min,Space_x_min), min(Cluster_x_max,Space_x_max),[],options);
end

CurrentEvals = CurrentEvals + FoodSize;

%% Choose Minimum Between Best Eagle and Best Food
[BestValue, FinalIndex] = sort([FitFood;NewFitFood; Fitness]);
% x = [NewEagle_Star;Food_Star;New_Food_Star];
% Star = x(BestIndex(1),:);

if NewFitFood < FitFood
    Food_Star = New_Food_Star;
    FitFood = NewFitFood;
% else
%     NewEagle_Star = Eagle_Star;
%     Fitness = bestold;
end

%% Plot Best Eagle/Food
% delete(ES); delete(S);
% hold on; ES = plot(NewEagle_Star(1), NewEagle_Star(2),'h','markersize',13,'LineWidth',1.5);
% hold on; S = plot(Star(1), Star(2),'*','markersize',13,'LineWidth',1.5); pause(0.1);

%% Record Best Eagle and Initial Best Food
bestcost = BestValue(1);
Result = [Result; G NewEagle_Star Fitness Food_Star FitFood bestcost CurrentEvals];

%% Check Stopping Criterion
% if CurrentEvals >= MaxEvals
%     break
% end

%% Update
Eagle_Star = NewEagle_Star;
bestold = Fitness;

end

%fprintf('Result: Gen # | Best Eagle x | Best Eagle y | Fitness | Best Food x | Best Food y | Fitness | Best Fitness | Func Evals \n');
%disp(Result)

toc

fbest_pheagle = Result(end,end-1);
if FinalIndex(1) == 3
    xbest_pheagle = NewEagle_Star;
else
    xbest_pheagle = Food_Star;
end
evalnum_pheagle = CurrentEvals;

end

% ___________________________________
function o=Levy(d)
betaL=1.5;
sigma=(gamma(1+betaL)*sin(pi*betaL/2)/(gamma((1+betaL)/2)*betaL*2^((betaL-1)/2)))^(1/betaL);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/betaL);
o=step;
end