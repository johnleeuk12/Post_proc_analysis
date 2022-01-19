function figures_tones(out)


N_list = out.data{1,1}(:,1);

animal_name = 'M60F';

% B : best dB, BF, peak lat, FR, min lat,


B = zeros(length(N_list),7);

for n = 1:length(N_list)
    FR = zeros(1,4);
    db_count = 0;
    for db_ind = 1:4
        if ~isempty(out.data{1,db_ind+1}{n,3})
            db_count = db_count +1;
            FR(db_ind) = out.data{1,db_ind+1}{n,3};
        else
            FR(db_ind) = 0;
        end
    end
    [~, bestdb_ind] = max(FR);
    if length(out.data{1,bestdb_ind+1}{n,1}) == 1
        if ~isempty(out.data{1,bestdb_ind+1}{n,5})
            B(n,:) = [bestdb_ind, out.data{1,bestdb_ind+1}{n,1}, out.data{1,bestdb_ind+1}{n,2},...
                out.data{1,bestdb_ind+1}{n,3},out.data{1,bestdb_ind+1}{n,4},out.data{1,bestdb_ind+1}{n,5},db_count];
        else
            B(n,:) = [bestdb_ind, out.data{1,bestdb_ind+1}{n,1}, out.data{1,bestdb_ind+1}{n,2},...
                out.data{1,bestdb_ind+1}{n,3},out.data{1,bestdb_ind+1}{n,4},NaN,db_count];
        end
    end
    if B(n,1) == 0
        B(n,2:6) = NaN;
        B(n,7) = 0;
    end
end


%%
global C
C = {};
C.pool.neurons_info = out.data{1,1};
C.pool.data = B;
C.pool.data_tags = {'db_ind'...
                    'BF, kHz'...
                    'peak_latency, ms'...
                    'peak_FR, spks/s'...
                    'min_latency, ms'};

H_list = unique(C.pool.neurons_info(:,2));
for h = 1:length(H_list)
    hind = C.pool.neurons_info(:,2) == H_list(h);
    htracks = unique(C.pool.neurons_info(hind,3));
    for ht = 1:length(htracks)
        %remove units without min_lat
        
        ht_ind = (C.pool.neurons_info(:,2) == H_list(h) & C.pool.neurons_info(:,3) == htracks(ht) &~isnan(B(:,6))); % & B(:,7)>2);
        ht_ind2 = (C.pool.neurons_info(:,2) == H_list(h) & C.pool.neurons_info(:,3) == htracks(ht));
        
        C.H{H_list(h)}{htracks(ht)}.BF.mean = mean(B(ht_ind,2),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.BF.med = median(B(ht_ind,2),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.BF.std = std(B(ht_ind,2),'omitnan');
        
        C.H{H_list(h)}{htracks(ht)}.pklat.mean = mean(B(ht_ind,3),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.pklat.med = median(B(ht_ind,3),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.pklat.std = std(B(ht_ind,3),'omitnan');
        
        C.H{H_list(h)}{htracks(ht)}.pkFR.mean = mean(B(ht_ind,4),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.pkFR.med = median(B(ht_ind,4),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.pkFR.std = std(B(ht_ind,4),'omitnan');
        
        C.H{H_list(h)}{htracks(ht)}.spontFR.mean = mean(B(ht_ind2,5),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.spontFR.med = median(B(ht_ind2,5),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.spontFR.std = std(B(ht_ind2,5),'omitnan');
        
        C.H{H_list(h)}{htracks(ht)}.minlat.mean = mean(B(ht_ind,6),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.minlat.med = median(B(ht_ind,6),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.minlat.std = std(B(ht_ind,6),'omitnan');
        
        C.H{H_list(h)}{htracks(ht)}.bestdB.mean = mean(B(ht_ind,7),'omitnan');
        C.H{H_list(h)}{htracks(ht)}.bestdB.median = median(B(ht_ind,7),'omitnan');

    end


end

%% PT analysis 

filedir = fullfile('E:\DATA',filesep,animal_name,filesep,'ana_tones\data');

load(fullfile(filedir,filesep,'neurons_loc_tag.mat'));

%%
N = length(neurons_loc_tag);
Core = {};
Belt = {};

Core.params.N = 0;
Core.params.out_ind = [];
Core.params.out_ind2 = [];
Core.params.out_ind3 = [];
Core.params.list = [5,6,15,16];
% Core.params.list = [2,3,6];


Belt.params.N = 0;
Belt.params.out_ind = [];
Belt.params.out_ind2 = [];
Belt.params.list = [4,10,11,14];
% Belt.params.list = [7:10];
% out_ind1 is for units with peak detected 
% out_ind2 is for units with min_lat. these two are different ways of
% detecting significant responses. 

% calculating number of PT responsive units
for h = Core.params.list
    Core.params.N = Core.params.N + length(find([neurons_loc_tag.hole_nb] == h));
    Core.params.out_ind = [Core.params.out_ind; find(C.pool.neurons_info(:,2) == h)];
    Core.params.out_ind2 = [Core.params.out_ind2; find(C.pool.neurons_info(:,2) == h & ~isnan(C.pool.data(:,6)))];
%     Core.params.out_ind2 = [Core.params.out_ind2; find(C.pool.neurons_info(:,2) == h & ~isnan(C.pool.data(:,6)) & C.pool.data(:,1)>2)];
end
% Core.params.out_ind2 = Core.params.out_ind;

for h = Belt.params.list
    Belt.params.N = Belt.params.N + length(find([neurons_loc_tag.hole_nb] == h));
    Belt.params.out_ind = [Belt.params.out_ind; find(C.pool.neurons_info(:,2) == h)];
    Belt.params.out_ind2 = [Belt.params.out_ind2; find(C.pool.neurons_info(:,2) == h & ~isnan(C.pool.data(:,6)))];
%     Belt.params.out_ind2 = [Belt.params.out_ind2; find(C.pool.neurons_info(:,2) == h & ~isnan(C.pool.data(:,6))& C.pool.data(:,1)>2)];
end

% Belt.params.out_ind2 = Belt.params.out_ind;

Belt.N_PT_resp =length(Belt.params.out_ind)-sum(C.pool.data(Belt.params.out_ind,1) ==0);
Core.N_PT_resp =length(Core.params.out_ind)-sum(C.pool.data(Core.params.out_ind,1) ==0);
Belt.N_PT_resp2 =length(Belt.params.out_ind2);
Core.N_PT_resp2 =length(Core.params.out_ind2);


%% distribution

%change k_ind to plot different properties

% edges = 0:0.1:10; %ms
% k_ind = 5;
for h = Core.params.list
figure
ind2 = find(out.data{1,1}(Core.params.out_ind2,2) ==h);
ind2= Core.params.out_ind2(ind2);
% subplot(1,2,1)
histogram(C.pool.data(ind2,k_ind)-200,edges)
% subplot(1,2,2)
% boxplot(C.pool.data(ind2,k_ind)-200);
% axis([0.5 1.5, 0 400])
sgtitle(['H', num2str(h),' median = ', num2str(median(C.pool.data(ind2,k_ind))-200),...
    'SD = ' num2str(std(C.pool.data(ind2,k_ind)))])
end

% figure
% boxplot(C.pool.data(Core.params.out_ind2,k_ind)-200,out.data{1,1}(Core.params.out_ind2,2))




%%

for h = Belt.params.list
figure
ind2 = find(out.data{1,1}(Belt.params.out_ind2,2) ==h);
ind2= Belt.params.out_ind2(ind2);
% subplot(1,2,1)
histogram(C.pool.data(ind2,k_ind)-200,edges)
% subplot(1,2,2)
% boxplot(C.pool.data(ind2,k_ind)-200);
% axis([0.5 1.5, 0 400])
sgtitle(['H', num2str(h),' median = ', num2str(median(C.pool.data(ind2,k_ind))-200),...
    'SD = ' num2str(std(C.pool.data(ind2,k_ind)))])
end


%%

edges = 0:10:200; %ms
k_ind =7;
figure(10)
histogram(C.pool.data(Core.params.out_ind2,k_ind),edges)
hold on
histogram(C.pool.data(Belt.params.out_ind2,k_ind),edges)

hold off
title([num2str(median(C.pool.data(Core.params.out_ind2,k_ind),'omitnan')),'      ', ...
    num2str(median(C.pool.data(Belt.params.out_ind2,k_ind),'omitnan'))])
% median(C.pool.data(Belt.params.out_ind2,k_ind))-200
% % end


disp(['mean A1:', num2str(mean(C.pool.data(Core.params.out_ind2,k_ind),'omitnan')),...
        '+', num2str(std(C.pool.data(Core.params.out_ind2,k_ind),'omitnan'))]);
disp(['mean A1:', num2str(mean(C.pool.data(Belt.params.out_ind2,k_ind),'omitnan')),...
        '+', num2str(std(C.pool.data(Belt.params.out_ind2,k_ind),'omitnan'))]);
    
p = ranksum(C.pool.data(Core.params.out_ind2,k_ind),C.pool.data(Belt.params.out_ind2,k_ind))





% 
figure(11)
boxplot(C.pool.data([Core.params.out_ind2; Belt.params.out_ind2],k_ind), ...
        [out.data{1,1}(Core.params.out_ind2,2); out.data{1,1}(Belt.params.out_ind2,2)+100])













































%%
% 
% X = {};
% 
% db_ind = 3; % 4 is 40db, 3 is 60db
% for n = 1:length(out.data{db_ind})
%     if ~isempty(out.data{db_ind}{n})
%         X(n).cf = out.data{db_ind}{n}(1);
%     else
%         X(n).cf = [];
%     end
%     X(n).Hole = out.data{6}(n,1);
%     X(n).Track = out.data{6}(n,2);
% end
% 
% %%
% 
% tracks = unique([X.Track]);
% edges = logspace(3,4.5,50);
% figure
% cf_mean = [];
% for t = 1:length(tracks)
%     ind = find([X.Track] == tracks(t));
%     cf_pool = [X(ind).cf];
%     histogram(cf_pool,edges)
%     hold on
%     cf_mean = [cf_mean;[mean(cf_pool),tracks(t)]];
% end
% 
% xticks([1e3 : 2.5*1e3:3*1e4])
% % set(gca,'xticklabel',num2str(get(gca,'xtick')/1e3','%.1f'))
% xticklabels(num2str(get(gca,'xtick')'/1e3))
% xlabel('kHz')
% set(gca,'xscale','log')