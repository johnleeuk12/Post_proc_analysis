
X = {};

db_ind = 3; % 4 is 40db, 3 is 60db
for n = 1:length(out.data{db_ind})
    if ~isempty(out.data{db_ind}{n})
        X(n).cf = out.data{db_ind}{n}(1);
    else
        X(n).cf = [];
    end
    X(n).Hole = out.data{6}(n,1);
    X(n).Track = out.data{6}(n,2);
end

%%

tracks = unique([X.Track]);
edges = logspace(3,5,100);
figure
cf_mean = [];
for t = 1:length(tracks)
    ind = find([X.Track] == tracks(t));
    cf_pool = [X(ind).cf];
    histogram(cf_pool,edges)
    hold on
    cf_mean = [cf_mean;[mean(cf_pool),tracks(t)]];
end

xticks([1e3 : 2.5*1e3:3*1e4])
% set(gca,'xticklabel',num2str(get(gca,'xtick')/1e3','%.1f'))
xticklabels(num2str(get(gca,'xtick')'/1e3))
xlabel('kHz')
set(gca,'xscale','log')
%%
% save_dir = 'E:\DATA\ana_tones\M60F_map';
% hole_number = X(1).Hole;
% if db_ind == 4
%     
%     save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_40dB.mat']),'X');
% elseif db_ind == 3
%     save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_60dB.mat']),'X');
% elseif db_ind == 2
%     
%     save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_80dB.mat']),'X');
% end
% 
% disp('file saved');

%% tuning PSTH per track

% stimrange = [2.^[-3:0.05:8]].';
% stimrange = stimrange*1e3; %from 500Hz to 32kHz, 20 steps/ octave


tracks = unique([X.Track]);
X_PSTH = {};
n_list_tracks = {};
track_cf = {};
ST = {};
figure
for t = 1:length(tracks)
    
    %     X_PSTH{t} = zeros(length(stimrange),600);
    X_PSTH{t} = [];
    n_list_tracks{t} = find(out.data{6}(:,2) == tracks(t));
    i = 1;
    stim_list = [];
    ST{t} = stim_list;
    try
    while isempty(stim_list)
        stim_list = out.data{1,db_ind}{n_list_tracks{t}(i),3};
        i = i+1;
    end
    ST{t} = stim_list;
    for n = 1:length(n_list_tracks{t})
        if length(out.data{1,db_ind}{n_list_tracks{t}(n),3} )==61 || length(out.data{1,db_ind}{n_list_tracks{t}(n),3} )==81 || length(out.data{1,db_ind}{n_list_tracks{t}(n),3} )==101
            %             first_Hz = out.data{1,db_ind}{n_list_tracks{t}(n),3}(1);
            %             last_Hz = out.data{1,db_ind}{n_list_tracks{t}(n),3}(end);
            if ~isempty(out.data{1,db_ind}{n_list_tracks{t}(n),2})
                if isempty(X_PSTH{t})
                    X_PSTH{t} = out.data{1,db_ind}{n_list_tracks{t}(n),2}(:,1:600);
                else
                    try
                    X_PSTH{t} = X_PSTH{t} + out.data{1,db_ind}{n_list_tracks{t}(n),2}(:,1:600);
                    catch
                    end
                end
            end
%             X3 = imgaussfilt(out.data{1,db_ind}{n_list_tracks{t}(n),2},[3,20],'Padding',0);
%             X3 = X3-mean2(X3(:,1:150));
%             imagesc(X3)
%             drawnow
%             pause
        end
    end
    
    if isempty(X_PSTH{t}) %no 61 or 81
        for n= 1:length(n_list_tracks{t})
            if length(out.data{1,db_ind}{n_list_tracks{t}(n),3} )==41
                %             first_Hz = out.data{1,db_ind}{n_list_tracks{t}(n),3}(1);
                %             last_Hz = out.data{1,db_ind}{n_list_tracks{t}(n),3}(end);
                if isempty(X_PSTH{t})
                    X_PSTH{t} = out.data{1,db_ind}{n_list_tracks{t}(n),2}(:,1:600);
                else
                    X_PSTH{t} = X_PSTH{t} + out.data{1,db_ind}{n_list_tracks{t}(n),2}(:,1:600);
                end
                %             X3 = imgaussfilt(out.data{1,db_ind}{n_list_tracks{t}(n),2},[3,20],'Padding',0);
                %             X3 = X3-mean2(X3(:,1:150));
                %             imagesc(X3)
                %             drawnow
                %             pause
            end
        end
    end
        
    
    X2 = imgaussfilt(X_PSTH{t},[3,10],'Padding',0); %,'Padding','circular');
    X2 = X2-mean2(X2(:,1:150)); % removing spont rate
    subplot(2,ceil(length(tracks)/2),t)
    imagesc(X2)
    drawnow()
    [BF_ind, pos, ~]  = find(X2 == max(max(X2(:,200:end))));
    
    % checking whether peak is significant
    %                     if mean2(X(:,200:end)) < 1 % average firing rate after onset is below 1
    %                         BF_ind = [];
    %                     if mean2(abs(X(:,200:end)-mean2(X(:,1:150))))<1
    if X2(BF_ind,pos) < std2(X2(:,1:150))*2
        BF_ind = [];
    end
    if ~isempty(BF_ind)
        track_cf(t).cf = stim_list(BF_ind);
        track_cf(t).FR = X2(BF_ind,pos);
    else
        track_cf(t).cf = [];
    end
    track_cf(t).tracks = tracks(t);
    
    catch
    end

    
end



% tt =2;
% X2 = imgaussfilt(X_PSTH{tt},[3,5],'Padding',0); %,'Padding','circular');
% X2 = X2-mean2(X2(:,1:150)); % removing spont rate
% figure
% imagesc(X2)


%%
D ={};
D.unit_cf = X;
D.track_cf = track_cf;
D.track_PSTH(:).PSTH = X_PSTH.';
D.track_PSTH(:).stim = ST.';
% D.track_PSTH = {D.track_PSTH, ST.'};
save_dir = 'E:\DATA\ana_tones\M60F_tracks';
hole_number = X(1).Hole;
if db_ind == 4
    
    save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_40dB.mat']),'D');
elseif db_ind == 3
    save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_60dB.mat']),'D');
elseif db_ind == 2
    
    save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_80dB.mat']),'D');
end

disp('file saved');




%% separate code for H3
% 
% X = {};
% for n = 1:length(Y)
%     if isnan(Y(n,1))
%         X(n).cf = [];
%     else
%         X(n).cf = Y(n,1);
%     end
%     X(n).Hole = Y(n,3);
%     X(n).Track = Y(n,4);
% end
% hole_number = 3;
%     save(fullfile(save_dir,['M60F_H',num2str(hole_number),'_40dB.mat']),'X');
