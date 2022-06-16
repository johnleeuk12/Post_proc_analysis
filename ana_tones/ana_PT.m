function ana_PT()

%% calculating number of units per track
% we can use this to study effect of impedance.
animal_name = 'M160E';

savedir = fullfile('E:\DATA',filesep,animal_name,filesep,'ana_tones\data');
addpath(savedir);
load([savedir '\neurons_loc_tag.mat']);
H_list = unique([neurons_loc_tag.hole_nb]).';

track_data = [];
T_list_H = [];
for ho = 1:length(H_list)
    tr_list = unique([neurons_loc_tag(find([neurons_loc_tag.hole_nb] == H_list(ho))).track_nb]);
    T_list_H = [T_list_H; length(tr_list)];
    for tr = tr_list
        track_data = [track_data; H_list(ho), tr, ...
                    length(find([neurons_loc_tag.hole_nb] == H_list(ho) & [neurons_loc_tag.track_nb] == tr))];
    end
end
edges = 5:10:100;
figure(1)
histogram(track_data(:,3),edges)
title(['median: ' num2str(median(track_data(:,3)))])


%% unit SNR distribution


animal_name = 'M56E';


load([animal_name,'_unit_list_new.mat']);
load([animal_name,'_neurons_list.mat']);


SNR_data = zeros(size(neurons_list,1),1);
for n = 1:size(neurons_list,1)
    SNR_bin = [];
    for bin = 1:size(neurons_list,2)
        if ~isnan(neurons_list(n,bin))
            if neurons_list(n,bin) ~=0
            SNR_bin = [SNR_bin, unit_list.data(neurons_list(n,bin),2)];
            end
        end
    end
    SNR_data(n) = mean(cell2mat(SNR_bin));
end
            
           
edges = 5:2:50;
figure(2)
histogram(SNR_data,edges)
title(['median: ' num2str(median(SNR_data))])
%% brain barrier and noisefloor

addpath('C:\Users\John.Lee\Documents\GitHub\analysis-tools');
filedir = 'E:\DATA\Other\2021-12-24_14-50-13';

% chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
%     22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
%     64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
%     38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 59 61 57 55 51 53 49 48 54 52 56 58 62 60 64 ...
    37 41 39 46 44 34 50 36 47 33 35 43 45 40 42 38];
data = {};
parfor ch = 1:64
    
    [data{ch},~,~] = load_open_ephys_data_faster([filedir filesep '100_CH' num2str(ch) '.continuous' ]);
end



%%
rawD = zeros(64,length(data{1}));
for ch = 1:length(chanMap)
    rawD(ch,:) = data{chanMap(ch)};
end


filtData = single([]);
[b,a] = butter(4, [0.0244 0.6104],'bandpass');
% [b,a] = butter(4, 500/fs,'low');
Nbuff = 32*32;
Nbatch = length(rawD(1,:))/Nbuff;
for batch = 1:Nbatch
    ipoint = (batch-1)*Nbuff +1;
    datr = filter(b,a,single(rawD(:,ipoint:ipoint+Nbuff-1)).');
    datr = flipud(datr);
    datr = filter(b,a,datr);
    datr = flipud(datr);
    filtData(:,ipoint:ipoint+Nbuff-1) = int16(datr.');
end

filtData = filtData-median(filtData(1:8,:));



filtData = double(filtData);



figure(4)
hold off
% rectangle('Position',[0 -7000 50000 7000],'FaceColor',[0,0,0],'EdgeColor','none')
hold on
for ch = 1:64
    plot(filtData(ch,1:50000) + ch*-75,'k');
    
end
% bad channels 44 47 55
            
stdD = std(filtData,0,2);
figure(5)
scatter([1:64],stdD)
% plot(stdD);
% meanD = mean(abs(filtData),2);
% plot(meanD);
%             
            
            
            
        %% 02/25/2022 Find multipeak responses
% db_ind.      4 is 20db % 3 is 40db, 2 is 60db 1 is 80 db


for n = 1:50
    
    fi = figure(n);
    set(fi, 'Position', [100 200 2400 660]);
    for db_ind = 2:5
        out.data{1,db_ind}{n,8} = {};
        if ~isempty(out.data{1,db_ind}{n,1})
            X = out.data{1,db_ind}{n,6};
            X2 = imgaussfilt(X,[3,20],'Padding',0); %,'Padding','circular');
            X2 = X2-mean2(X2(:,1:150)); % removing spont rate
            subplot(1,4,db_ind-1)
            imagesc(X2)
            hold on
            rectangle('Position',[200 0 100 size(X2,1)],'LineWidth',2)
            p = 1;
            
            for tf = 201:100:301
                T = mean(X2(:,tf:tf+100),2);
                SD = std2(X2(:,1:150));
                [pks, locs] = findpeaks(T,'MinPeakProminence',2*SD);
                locs2 = [];
                if ~isempty(locs)
                    for lc = 1:length(locs)
                        if pks(lc)>4*SD
                            scatter(tf+50,locs(lc),100,'rp','filled')
                            locs2 = [locs2 locs(lc)];
                        end
                    end
                end
                out.data{1,db_ind}{n,8}{p} = locs2;
                 p = p+1;
            end

        end
    end
    
    drawnow()
    pause(0.1)
    
    
    
end





    
            
            
            
            
            
            
            
