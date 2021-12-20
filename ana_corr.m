% 12/2/2021 correlation analysis
%{
Here we attempt to calculate the correlation of time-variant activity
(during stim and during offset) between a given units firing patterns across
stimuli
%}

mean_thresh= 1;
m = mean([DD.data],2)*1e3;
keep_neurons = m >= mean_thresh;
nreps = 8;
nid_list = good_list_nid(keep_neurons);
TD = {};

%% Trial average D 

%%
% unsmoothened TD %% Wrong data! average response correlation, not
% individual unit correlation
% for st = 1:stim_num
%     TD(st).data = zeros(sum(keep_neurons),PreStim+PostStim+stim_dur);
%     for nr = 1:nreps
%         TD(st).data = TD(st).data +mean(DD((st-1)*nreps+nr).data(keep_neurons,:));
%     end
% end

% %smoothened TD
binwidth = 20; %ms
PreStim = floor(PreStim/binwidth);
PostStim = floor(PostStim/binwidth);
stim_dur = floor(stim_dur/binwidth);

%%
smooth_bin = 20;
% if smooth_bin ~= 1
%     for i = 1:length(DD)
%         DD(i).data = smoother(DD(i).data, smooth_bin , 1);
%     end
% end


for st = 1:stim_num
    TD(st).data = zeros(sum(keep_neurons),PreStim+PostStim+stim_dur);
    for nr = 1:nreps
        for b = 1:PreStim+PostStim+stim_dur
            TD(st).data(:,b) = TD(st).data(:,b) +mean(DD((st-1)*nreps+nr).data(keep_neurons,(b-1)*binwidth+1:b*binwidth),2);
        end
    end
    TD(st).data = smoother(TD(st).data,smooth_bin,binwidth);
end
%%
% figure(100)
% for n = 1:200
%     for st = 1:stim_num
%     plot(TD(st).data(n,:)*1e3);
%     hold on
%     end
%     pause
%     hold off
% end


%% Correlation matrix onset or offset
tpoint = [PreStim,PreStim+stim_dur];
% tpoint = [0,PreStim];
% tpoint = [PreStim+stim_dur,PreStim+stim_dur+PostStim];
CorrMatrix = zeros(stim_num);
new_nid_count = 0;
new_nid_list = [];
for n = 1:length(nid_list)
    coorD = zeros(stim_num,tpoint(2)-tpoint(1));
    remove_n = 0;
    for st = 1:stim_num
        coorD(st,:) = TD(st).data(n,tpoint(1)+1:tpoint(2));
    end
    mean_D = mean(coorD,2);
    for st = 1:stim_num
        if mean_D(st)==0
            remove_n = 1;
        end
    end
    % average response during stim
    
    if remove_n==0
        [R,P] = corrcoef(coorD.','Rows','complete');
        %           [R,P] = corrcoef(mean_D.','Rows','complete');
        for r = 1:length(R)^2
            if P(r) > 0.05   && R(r) ~= 1
                R(r) = 0;
            end
        end
        CorrMatrix = CorrMatrix + R;
        new_nid_count = new_nid_count+1;
        new_nid_list = [new_nid_list nid_list(n)];
        PoolCorr(:,:,new_nid_count) = R;
        
        
    end
end

figure(56)
CorrMatrix = CorrMatrix/new_nid_count;
colormap(flipud(hot))
% imagesc(stim_set,[-5:0.25:5],CorrMatrix)
imagesc(stim_set,[-5:0.25:5],CorrMatrix);
caxis([0 0.3])

%%
% mean correlation per stim

mean_corr = zeros(1,stim_num);
error_corr = zeros(1,stim_num);

for st = 1:stim_num
%     mean_corr(st) = (sum(CorrMatrix(:,st))-1)/(stim_num-1);
mean_corr(st) = median(CorrMatrix(:,st));
    if st == 1
        error_corr(st) = std(CorrMatrix(st+1:end,st),0,1);
    elseif st == stim_num
        error_corr(st) = std(CorrMatrix(1:st-1,st),0,1);
    else
         error_corr(st) = std([CorrMatrix(1:st-1,st);CorrMatrix(st+1:end,st)],0,1);
    end
end

error_corr1 = error_corr/stim_num;

figure(58)
% errorbar([-1:0.5:5], mean_corr,error_corr1)
errorbar([-5:0.25:5],movmean(mean_corr,3),movmean(error_corr1,3));
hold on
% axis([-5 1 0.1 0.15])

%% average correlation for n, n-1 and n+1 SD
% 
for r = 1:length(CorrMatrix)^2
    if CorrMatrix(r) == 1
        CorrMatrix(r) = NaN;
    end
end
%         CorrMatrix(r) = (CorrMatrix(r-1)+CorrMatrix(r+1)/2

% % for trills
% mean_corr_group = zeros(5,1);
% error_corr_group = zeros(5,1);
% for sd = 0:4
%     sub_corr = CorrMatrix(sd*2-2+3:sd*2+2+3,sd*2-2+3:sd*2+2+3);
%     mean_corr_group(sd+1) =  mean(mean(sub_corr,'omitnan'));
%     error_corr_group(sd+1) = std(reshape(sub_corr,[],1),'omitnan')/length(sub_corr);
% end

%for phees
    
mean_corr_group = zeros(9,1);
error_corr_group= zeros(9,1);
sdd = 1;
for sd = 5:4:37
    
    sub_corr = CorrMatrix(sd-4:sd+4,sd-4:sd+4);
    mean_corr_group(sdd) = mean(mean(sub_corr,'omitnan'));
    error_corr_group(sdd) = std(reshape(sub_corr,[],1),'omitnan')/length(sub_corr);
    sdd =sdd+1;
end

figure(59)
% subplot(1,2,2)
% errorbar([-4:4],mean_corr_group/max(mean_corr_group),error_corr_group/max(mean_corr_group))
% hold on
% subplot(1,2,1)
errorbar([-4:4],mean_corr_group,error_corr_group)
hold on
% plot(mean_corr_group)



%% find units with high local correlation 

c = hot(9);
sdd = 0;
for sd = 5:4:21 %37
    sub_corr_mean = mean(mean(PoolCorr(sd-4:sd+4,sd-4:sd+4,:),1,'omitnan'),2);
    sub_corr_mean = reshape(sub_corr_mean,1,[]);
    [S, I] = sort(sub_corr_mean);
    sdd = sdd+1;
    figure(100)
    hold on
    plot(S,'LineWidth',2,'color',c(sdd,:))
end
hold off


%% average corr per SD
% stim_ind = 1:41;
% g = ceil(abs(stim_ind-21)/2);
% 
% mean_corr2 = zeros(1,max(g));
% error_corr2 = zeros(1,max(g));
% for st = 1:max(g)
%     ind_g = find(g ==st);
%     mean_corr2(st) = mean(mean_corr(ind_g));
%     error_corr2(st) = mean(error_corr(ind_g))/sqrt(stim_num*length(ind_g));
% end
%     
% 
% 
% figure(57)
% errorbar([0.5:0.5:5],mean_corr2,error_corr2);
% axis([0.5 5 0.048 0.065])
% 
% 
% 
% 
% 
% %% Trills
% 
% 
% figure(101)
% for st = 1:stim_num
%     scatter([-1:0.5:5],CorrMatrix(st,:));
%     hold on 
% end
% 
% 



















