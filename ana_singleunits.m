%% find SU index for selected nids
nid_list = good_list_nid(keep_neurons);
Sel_nid = nid_list; % try all units

NN = length(SUrate);
SU_ind= zeros(1,length(Sel_nid));
i =0;
for n = Sel_nid
    i = i+1;
    for nn = 1:NN
        if SUrate{nn}{1}.nid == n
            SU_ind(i) = nn;
        end
    end
end

switch VT_code
    case {10, 11}
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,10);
        end
    case 01
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,8);
        end
end
%% 2. select units with significant response 

tpoint3 = PreStim + stim_dur+ 100;
tpoint4 = tpoint3 + 200;
% tpoint1 = PreStim;
% tpoint2 = tpoint1 + 200;
tpoint1 = PreStim+100;
tpoint2 = tpoint1 + stim_dur;
tdur1 = tpoint2-tpoint1;
tdur2 = tpoint4-tpoint3;
mat_tuning = [];
mat_tuning2 = [];
mat_peak_ind = [];
Sel_nid_all = [];
DR_peak = [];
spont_thresh = [];


for n = 1:length(SU_ind)
    
% % %     for trills
% %     DR = zeros(length(unique(SUrate{SU_ind(n)}{1}.stim(:,1))),tdur+1);
% %     
% %     stim_ind = find(SUrate{SU_ind(n)}{1}.stim(:,2) < 30);
% %     stim_set = SUrate{SU_ind(n)}{1}.stim(stim_ind,1);
% %     for st = 1:length(stim_ind)
% %         DR(st,:) = mean(SUrate{SU_ind(n)}{1}.PSTH{st}(:,tpoint1:tpoint2));
% %     end
    DR = zeros(length(SUrate{SU_ind(n)}{1}.mean),tdur1+1);
    DR2 = zeros(length(SUrate{SU_ind(n)}{1}.mean),tdur2+1);
    
    for st = 1:length(SUrate{SU_ind(n)}{1}.mean)
        DR(st,:) = mean(SUrate{SU_ind(n)}{1}.PSTH{st}(:,tpoint1:tpoint2)); %change which is being sorted 
        DR2(st,:) = mean(SUrate{SU_ind(n)}{1}.PSTH{st}(:,tpoint3:tpoint4)); 
    end
    
    DR_mean = mean(DR,2);
%     DR_mean = movmean(DR_mean,3);
    [~, ind] = max(DR_mean);
    std_thresh = mean(SUrate{SU_ind(n)}{1}.spont) + 2*std(SUrate{SU_ind(n)}{1}.spont);
    if ind == length(DR_mean)
        ind = ind-1;
    elseif ind == 1
        ind = ind +1; 
    end
       
    
    if DR_mean(ind-1:ind+1) > std_thresh   
%     if DR_mean(ind)>std_thresh
        mat_tuning = [mat_tuning;(mean(DR,2)-mean(SUrate{SU_ind(n)}{1}.spont)).'];
        mat_tuning2 = [mat_tuning2;(mean(DR2,2)-mean(SUrate{SU_ind(n)}{1}.spont)).'];
        [DP, real_ind] = max(DR_mean);
        DR_peak = [DR_peak;[DP, (DP-DR_mean(21))/range(DR_mean)]];
        mat_peak_ind = [mat_peak_ind; real_ind];
        Sel_nid_all = [Sel_nid_all; SUrate{SU_ind(n)}{1}.nid];
        spont_thresh = [spont_thresh; std_thresh];
    end
end

peak_ind2 = abs(mat_peak_ind-20);
% scatter(peak_ind2,DR_peak(:,2))
% DR_peak3 = DR_peak(:,2)./peak_ind2;
% scatter(peak_ind2,DR_peak3)
% axis([-inf, inf, 0, 0.4])

%%
[B,I] = sort(mat_peak_ind);

sorted_mat_tuning = mat_tuning(I,:);
sorted_mat_tuning2 = mat_tuning2(I,:);
% soft normalization
mat2_peak_ind = [];
for n = 1:size(sorted_mat_tuning,1)
    [~, ind] = max(sorted_mat_tuning2(n,:));
    std_thresh = spont_thresh(I(n));
    
    
    if ind == length(DR_mean)
        ind = ind-1;
    elseif ind == 1
        ind = ind +1;
    end
    
    if sorted_mat_tuning2(n,ind-1:ind+1) > std_thresh
        [~, real_ind] = max(sorted_mat_tuning2(n,:));
        mat2_peak_ind = [mat2_peak_ind;real_ind];
    else
        mat2_peak_ind = [mat2_peak_ind;nan];
    
    end
    
end


for n = 1:size(sorted_mat_tuning,1)
    sorted_mat_tuning(n,:) = sorted_mat_tuning(n,:)/(range(sorted_mat_tuning(n,:))+5);
    sorted_mat_tuning2(n,:) = sorted_mat_tuning2(n,:)/(range(sorted_mat_tuning2(n,:))+5);
end


figure(58)
% mat_tuning2= imgaussfilt(sorted_mat_tuning,[0.5,1],'Padding',0);
% imagesc(SUrate{n}{1}.stim(:,1),[1:size(sorted_mat_tuning,1)],mat_tuning2)

% for phees
imagesc(SUrate{n}{1}.stim(:,1),[1:size(sorted_mat_tuning,1)],sorted_mat_tuning)
caxis([-1, 1]);
colormap(parula);
figure(61)
imagesc(SUrate{n}{1}.stim(:,1),[1:size(sorted_mat_tuning2,1)],sorted_mat_tuning2)
caxis([-1, 1]);
colormap(parula);


figure(59)
stimB = SUrate{SU_ind(n)}{1}.stim(B,1);
scatter(SUrate{SU_ind(n)}{1}.stim(B,1),[1:size(sorted_mat_tuning,1)]/size(sorted_mat_tuning,1),'.g')
% scatter(SUrate{SU_ind(n)}{1}.stim(B,1),[1:size(sorted_mat_tuning,1)],'.r')

hold on
nan_count = 0;
for sc = 1:length(mat2_peak_ind)
    if ~isnan(mat2_peak_ind(sc))
        scatter(SUrate{SU_ind(n)}{1}.stim(mat2_peak_ind(sc),1),sc/size(sorted_mat_tuning,1),'.r')
    
    else
        nan_count = nan_count +1;
    end    
end
% find peak for mat_tuning2
title([num2str(length(find(B>12 & B< 30))/length(B))])



% % for trills
% imagesc(stim_set,[1:size(sorted_mat_tuning,1)],sorted_mat_tuning)
% 
% figure(59)
% stimB = stim_set(B,1);
% scatter(stim_set(B,1),[1:size(sorted_mat_tuning,1)]/size(sorted_mat_tuning,1),'.b')
% hold on




%% compare PC and units
s_ind = zeros(length(Sel_nid_all),1);
for s = 1:length(Sel_nid_all)
    s_ind(s) =  find(nid_list == Sel_nid_all(s));
end

s_lat_C_value = lat_C(s_ind);
sum(lat_C(s_ind))/sum(lat_C);

%% save PSTH etc from nid

nid_list2 = intersect(nid_offset_phee,nid_offset_PT);

for n = nid_list2
    i = i+1;
    for nn = 1:NN
        if SUrate{nn}{1}.nid == n
            unit_data.PSTH = SUrate{nn}{1}.PSTH;
            unit_data.stim = SUrate{nn}{1}.stim;
        end
    end
end



%% calculate distance for tuning

tuning_peak = SUrate{SU_ind(n)}{1}.stim(mat_peak_ind,1);
tuning_peak = abs(tuning_peak-7.16);

figure(56)
edges = [0:0.12:2.4];
histogram(tuning_peak,edges)

%% response curves and magnet effect? 


g = ceil(peak_ind2/4);
for p = 1:length(g)
    if g(p) == 0
        g(p) = nan;
    elseif g(p) == 6
        g(p) = 5;
    end
end
figure(56)
boxplot(DR_peak(:,2),g)
axis([0.5,5.5,-0.05,1.05])



























% 
% 
% 
% % run after n.2
% 
% % NN = length(SUrate);
% % tpoint1 = PreStim + stim_dur+ 100;
% % tpoint2 = tpoint1 + 400 ;
% % tdur = tpoint2-tpoint1;
% 
% 
% for nn = Sel_nid_all.'
%     %     nn = Sel_nid(21);
%     
%     figure(50)
%     for n = 1:NN
%         if SUrate{n}{1}.nid == nn
%             
%             DR = zeros(length(SUrate{n}{1}.mean),tdur+1);
%             for st = 1:length(SUrate{n}{1}.mean)
%                 DR(st,:) = mean(SUrate{n}{1}.PSTH{st}(:,tpoint1:tpoint2));
%             end
%             
%             
%             nPSTH = zeros(length(SUrate{n}{1}.mean),trial_dur);
%             for st = 1:length(SUrate{n}{1}.mean)
%                 nPSTH(st,:) = mean(SUrate{n}{1}.PSTH{st}(:,1:trial_dur));
%             end
%             
%             PSTH2 = imgaussfilt(nPSTH,[3,20],'Padding',0);
%             subplot(1,2,1)
%             imagesc([1:trial_dur]-PreStim,SUrate{n}{1}.stim(:,1),PSTH2);
% 
%             colormap(flipud(gray))
%             xline(0);
%             xline(stim_dur);
%             ylabel('kHz');
%             xlabel('ms');
%             title(['M60Fu' num2str(nn)]);
%             subplot(1,2,2)
%             DR_2 = mean(DR,2)-mean(SUrate{n}{1}.spont);
%             errorbar(SUrate{n}{1}.stim,DR_2,std(DR,0,2)/sqrt(300))
% %             if DR_2(5)>1*std(SUrate{n}{1}.spont)
% %                 distance = abs(DR_2-DR_2(5));
% %                 subplot(2,2,3)
% %                 plot(stim_set,distance);
% %             end
% 
%         end
%     end
%     pause
%     clf
% end










































