% [D_phee, good_list_nid_phee,good_list_phee] = format4DataHigh(SUrate,10,1);


%% getting phee and PT 
keep_ind_phee = find(keep_neurons ==1);
nid_list_phee = nid_list;
aid_list_phee = aid_list;
D_phee = DD;

keep_ind_PT = find(keep_neurons ==1);
nid_list_PT = nid_list;
aid_list_PT = aid_list;
D_PT = DD;

%%
keep_ind_phee2 = [];
keep_ind_PT2 = [];
na_id_both = [];
for n = 1:length(nid_list_PT)
    ind = find(nid_list_phee == nid_list_PT(n) & aid_list_phee == aid_list_PT(n));
    na_id_both = [na_id_both, [nid_list_phee(ind);aid_list_phee(ind)]];
end

for n = 1:length(na_id_both)
    ind1 = find(nid_list_PT == na_id_both(1,n) & aid_list_PT == na_id_both(2,n));
    ind2 = find(nid_list_phee == na_id_both(1,n) & aid_list_phee == na_id_both(2,n));
    keep_ind_PT2 = [keep_ind_PT2; keep_ind_PT(ind1)];
    keep_ind_phee2 = [keep_ind_phee2; keep_ind_phee(ind2)];
end

for itrial = 1:length(D_phee)
    D_phee(itrial).data = D_phee(itrial).data(keep_ind_phee2,:);
end

for itrial = 1:length(D_PT)
    D_PT(itrial).data = D_PT(itrial).data(keep_ind_PT2,:);
end
%%
save('E:\DATA\Combined\ana_decoding\A1_PTvsPhee.mat','-v7.3','-nocompression');


%% normalizing 
ND = zeros(length(na_id_both),stim_num);
spont_phee = zeros(length(na_id_both),length(D_phee));
tpoints = [1,PreStim+20, PreStim+500];


for c = 1:length(D_phee)
    spont_phee(:,c) = mean(D_phee(c).data(:,tpoints(1):tpoints(2)),2);
end

ND_phee = zeros(stim_num,length(na_id_both),tpoints(3)-tpoints(2));
for st = 1:stim_num
    tempD = zeros(length(na_id_both),tpoints(3)-tpoints(2));
    for r = 1:nreps
        tempD = tempD + D_phee((st-1)*nreps+r).data(:,tpoints(2)+1:tpoints(3));
    end
    ND_phee(st,:,:) = tempD/nreps;
end

spont_phee = spont_phee*1e3;
ND_phee = ND_phee*1e3;

% ND = ND*1e3;
% for n = 1:size(ND,1)
%     ND(n,:) = ND(n,:)/(range(ND(n,:))+5);
% end
% ND_phee = ND.';


%% normalizing PT

% ND = zeros(length(na_id_both),length(D_PT));
spont_PT = zeros(length(na_id_both),length(D_PT));
% tpoints = [1,[D_PT(1).epochStarts(2:3)]];
tpoints = [1,PreStim+20, PreStim+500];

for c = 1:length(D_PT)
    spont_PT(:,c) = mean(D_PT(c).data(:,tpoints(1):tpoints(2)),2);
end

ND_PT = zeros(stim_num,length(na_id_both),tpoints(3)-tpoints(2));
for st = 1:stim_num
    tempD = zeros(length(na_id_both),tpoints(3)-tpoints(2));
    for r = 1:nreps
        tempD = tempD + D_PT((st-1)*nreps+r).data(:,tpoints(2)+1:tpoints(3));
    end
    ND_PT(st,:,:) = tempD/nreps;
end



spont_PT = spont_PT*1e3;
ND_PT = ND_PT*1e3;
% for n = 1:size(ND,1)
%     ND(n,:) = ND(n,:)/(range(ND(n,:))+5);
% end
% ND_PT = ND.';




%%

ND_phee_mean = mean(ND_phee,3)-mean(spont_phee,2).';
ND_PT_mean = mean(ND_PT,3)-mean(spont_PT,2).';

T = [];
for n = 1:length(na_id_both)
    [max_phee,max_ind_phee] = max(movmean(ND_phee_mean(:,n),3));
    [max_PT,max_ind_PT] = max(movmean(ND_PT_mean(:,n),3));
    if ND_phee_mean(max_ind_phee,n)> 1.5*std(spont_phee(n,:)) && ND_PT_mean(max_ind_PT,n)> 1.5*std(spont_PT(n,:))
        %         figure(n)
        %         plot(movmean(ND_phee_mean(:,n),3))
        %         hold on
        %         plot(movmean(ND_PT_mean(:,n),3))
        %         hold off
        %
                T = [T;n,max_ind_phee-max_ind_PT,max_ind_PT];
%         T = [T;n,(max_phee-max_PT)/max_PT,max_ind_PT];
    end
end

scatter(T(:,3),T(:,2))
























%% normalizing 
ND = zeros(length(na_id_both),length(D_phee));
spont = zeros(length(na_id_both),length(D_phee));
tpoints = [1,PreStim+20, PreStim+500];


for c = 1:length(D_phee)
    spont(:,c) = mean(D_phee(c).data(:,tpoints(1):tpoints(2)),2);
end

for c = 1:length(D_phee)
    ND(:,c) = mean(D_phee(c).data(:,tpoints(2):tpoints(3)),2)-mean(spont,2);
end



ND = ND*1e3;
% for n = 1:size(ND,1)
%     ND(n,:) = ND(n,:)/(range(ND(n,:))+5);
% end
ND_phee = ND.';


%% normalizing PT

ND = zeros(length(na_id_both),length(D_PT));
spont = zeros(length(na_id_both),length(D_PT));
% tpoints = [1,[D_PT(1).epochStarts(2:3)]];
tpoints = [1,PreStim+20, PreStim+500];

for c = 1:length(D_PT)
    spont(:,c) = mean(D_PT(c).data(:,tpoints(1):tpoints(2)),2);
end

for c = 1:length(D_PT)
    ND(:,c) = mean(D_PT(c).data(:,tpoints(2):tpoints(3)),2)-mean(spont,2);
end



ND = ND*1e3;
% for n = 1:size(ND,1)
%     ND(n,:) = ND(n,:)/(range(ND(n,:))+5);
% end
ND_PT = ND.';

%% Comparing responses during the first 500ms for PHees and PT. 








%% GLM training

n_bin = 5;
for s = 1:7
    stim_dur = 1180;
    stim_set = 4.76:0.12:9.56;
    stim_num = 41;
    zero_ind = 21;
    nreps = 10;
    %change here to choose stim pairs
    st_ind = 3+(s-1)*6;
    
    TD = zeros(length(na_id_both),nreps*n_bin*2);
    for r = 1:nreps*n_bin
        %         tempD = D((st_ind(st)-1)*nreps+r).data(:,tpoints(2):tpoints(3));
        %         TD(:,(st-1)*nreps+r) = mean(tempD,2);
        
        TD(:,r) = ND_PT(st_ind*nreps+r,:);
        TD(:,r+n_bin*nreps) = ND_phee(st_ind*nreps+r,:);
        
    end
    
    
    Y = [zeros(1,nreps*n_bin),ones(1,nreps*n_bin)];
    TD = TD.';
    Y = Y.';
    
    train_ind = [randsample(nreps*n_bin,7*n_bin); randsample(nreps*n_bin,7*n_bin)+10*n_bin];
    test_ind = setdiff(1:n_bin*2*nreps,train_ind);
    
    
    
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






































