% 12/15/2021 Across unit correlation analysis
%{
Here we attempt to calculate the correlation of time-variant activity
(during stim and during offset) across units firing patterns across
stimuli
%}

mean_thresh= 1;
m = mean([DD.data],2)*1e3;
keep_neurons = m >= mean_thresh;
nreps = 8;
nid_list = good_list_nid(keep_neurons);
TD = [];

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
    st_data = zeros(sum(keep_neurons),PreStim+PostStim+stim_dur);
    for nr = 1:nreps
        for b = 1:PreStim+PostStim+stim_dur
            st_data(:,b) = st_data(:,b) +mean(DD((st-1)*nreps+nr).data(keep_neurons,(b-1)*binwidth+1:b*binwidth),2);
        end
    end
    st_data = smoother(st_data,smooth_bin,binwidth);
    TD(st,:,:) = st_data;
end



%% average correlation


tpoint = [PreStim,PreStim+stim_dur];
% tpoint = [0,PreStim];
% tpoint = [PreStim+stim_dur,PreStim+stim_dur+PostStim];
TD_av = mean(TD(:,:,tpoint(1):tpoint(2)),3);

[CorrMatrix,P] = corrcoef(TD_av.','Rows','complete');
for r = 1:length(CorrMatrix)^2
    if P(r) > 0.05 && CorrMatrix(r) ~=1
        CorrMatrix(r) = 0;
    end
end

figure(56)
colormap(flipud(hot))
% imagesc(stim_set,[-5:0.25:5],CorrMatrix)
imagesc(stim_set,[-1:0.5:5],CorrMatrix);




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
errorbar([-1:0.5:5],movmean(mean_corr,3),movmean(error_corr1,3));
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

% for trills
mean_corr_group = zeros(5,1);
error_corr_group = zeros(5,1);
for sd = 0:4
    sub_corr = CorrMatrix(sd*2-2+3:sd*2+2+3,sd*2-2+3:sd*2+2+3);
    mean_corr_group(sd+1) =  mean(mean(sub_corr,'omitnan'));
    error_corr_group(sd+1) = std(reshape(sub_corr,[],1),'omitnan')/length(sub_corr);
end

%for phees
    
% mean_corr_group = zeros(9,1);
% error_corr_group= zeros(9,1);
% sdd = 1;
% for sd = 5:4:37
%     
%     sub_corr = CorrMatrix(sd-4:sd+4,sd-4:sd+4);
%     mean_corr_group(sdd) = mean(mean(sub_corr,'omitnan'));
%     error_corr_group(sdd) = std(reshape(sub_corr,[],1),'omitnan')/length(sub_corr);
%     sdd =sdd+1;
% end

figure(59)
% subplot(1,2,2)
% errorbar([-4:4],mean_corr_group/max(mean_corr_group),error_corr_group/max(mean_corr_group))
% hold on
% subplot(1,2,1)
errorbar([0:4],mean_corr_group,error_corr_group)
hold on
% plot(mean_corr_group)
% caxis([0 0.5])

%%
% %% temporal sequence correlation
% 
% % for st = 1:stim_num 
% st = 37;
% [R,P] = corrcoef(reshape(TD(st,:,tpoint(1):tpoint(2)),size(TD,2),[]).','Rows','complete');
% 
% % imagesc(R)
% for r = 1:length(R)^2
% %     if P(r)>0.05 && R(r) ~=1
% %         R(r) = 0;
% %     end
%     if isnan(R(r))
%         R(r) = 0;
%     end
% end
% cgs2 = clustergram(R);
% set(cgs2,'Linkage','complete','Dendrogram',15)
% % Y = pdist(R);
% 



%%
% CorrMatrix = zeros(stim_num,stim_num);
% for t = tpoint(1):tpoint(2)
%     [R,P] = corrcoef(TD(:,:,t).','Rows','complete');
%     for r = 1:length(R)^2
%         if P(r) > 0.05 && R(r) ~=1
%             R(r) = 0;
%         end
%     end
%     CorrMatrix = CorrMatrix + R;
% end
% CorrMatrix = CorrMatrix/size(TD,2);







