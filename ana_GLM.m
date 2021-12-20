%% GLM analysis 

D = DD;
% remove low firing rate neurons
mean_thresh= 1;
m = mean([D.data],2)*1e3;
keep_neurons = m >= mean_thresh;



for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(keep_neurons,:);
end

tpoints = [1,[D(1).epochStarts(2:3)]];


%% normalizing 
ND = zeros(sum(keep_neurons),length(D));
spont = zeros(sum(keep_neurons),length(D));

for c = 1:length(D)
    spont(:,c) = mean(D(c).data(:,tpoints(1):tpoints(2)),2);
end

for c = 1:length(D)
    ND(:,c) = mean(D(c).data(:,tpoints(2):tpoints(3)),2)-mean(spont,2);
end

ND = ND*1e3;
for n = 1:size(ND,1)
    ND(n,:) = ND(n,:)/(range(ND(n,:))+5);
end
ND = ND.';

% 
% figure
% for n = 1:size(ND,1)
%     plot(TD(n,:))
%     drawnow
%         pause
% 
% end

%%
B1 = {};
FitInfo1 = {};
yhat1 = {};
for s = 1:8
    stim_dur = 1180;
    stim_set = 4.76:0.12:9.56;
    stim_num = 41;
    zero_ind = 21;
    nreps = 10;
    %change here to choose stim pairs
    st_ind = [4+(s-1)*4,4+s*4];
    
    TD = zeros(sum(keep_neurons),2*nreps*3);
    for st = 1:2
        for r = 1:nreps*3
            %         tempD = D((st_ind(st)-1)*nreps+r).data(:,tpoints(2):tpoints(3));
            %         TD(:,(st-1)*nreps+r) = mean(tempD,2);
            
            TD(:,(st-1)*nreps*3+r) = ND((st_ind(st)-1)*nreps+r,:);
            
        end
    end
    
    % TD = TD*1e3;
    Y = [zeros(1,nreps*3),ones(1,nreps*3)];
    
    
    
    
    TD = TD.';
    Y = Y.';
    
    train_ind = [randsample(nreps*3,7*3); randsample(nreps*3,7*3)+10*3];
    test_ind = setdiff(1:60,train_ind);
    
    
    
    % X_train = TD(train_ind,:);
    % Y_train = Y(train_ind,:);
    X_test = TD(test_ind,:);
    
    % TD = [TD zeros(60,1)];
    
    
    
    [B1{s}, FitInfo1{s}] = lassoglm(TD, Y,'binomial','CV',10,'Alpha',1);
    B0 = FitInfo1{s}.Intercept;
    coef = [B0;B1{s}];
    yhat1{s} = glmval(coef,X_test,'logit');
    
    lassoPlot(B1{s},FitInfo1{s},'PlotType','CV');
    legend('show','Location','best') % show legend
    pause(0.1)
    drawnow
    disp(s)
end
%%


mincoefs = find(B1(:,idxLambdaMinDeviance))








lassoPlot(B1,FitInfo1,'PlotType','Lambda','XScale','log');
%%

figure(58)
p = 0;
for a = 0.1:0.1:1
    p = p+1;
    [B1, FitInfo1] = lassoglm(TD, Y,'binomial','NumLambda',25,'CV',10,'Alpha',a);
    subplot(2,5,p);
    lassoPlot(B1,FitInfo1,'PlotType','CV');
    legend('show','Location','best') % show legend
    title(['alpha= ' num2str(a)])
    drawnow
    pause(0.1)
end













