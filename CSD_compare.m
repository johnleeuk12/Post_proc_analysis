% Comparing CSD analysis and assigning spikes to depth/ layers
%12/02/2020


addpath('D:\DATA\M12E\Units'); % path to units
addpath('D:\DATA\M12E\Experiments'); %path to xbz files
load('M12E_neurons_list.mat');


if ~exist('M12E_unit_list_new.mat')
    load('M12E_unit_list.mat');
    fprintf('Time %3.0fs. adding analysis code... \n', toc);
    
    for i  = 1:length(M12E_unit_list.data)
        
        x = eval(M12E_unit_list.data{i,4});
        M12E_unit_list.data{i,6} = x.analysis_code;
        M12E_unit_list.data{i,7} = x.analysis_type;
        M12E_unit_list.data{i,8} = x.hole_number;
        M12E_unit_list.data{i,9} = x.track_number;
    end
    fprintf('Time %3.0fs. adding analysis code... Done! \n', toc);
    
else
    load('M12E_unit_list_new.mat');
end

%% The following code can be used universally, to extract spike and stim times from unit data.
% dB = 40;
% a_code = 230;
a_code = 1;
dB = 40;
hole_number = 3;
track_number = -1;
if size(hole_number) == 1
    u_list = find([M12E_unit_list.data{:,6}] == a_code & [M12E_unit_list.data{:,8}] == hole_number & [M12E_unit_list.data{:,9}] == track_number);
else
%     u_list = [];
%     for h = 1:length(hole_number)
%         u_list = [u_list find([M12E_unit_list.data{:,6}] == a_code & [M12E_unit_list.data{:,8}] == hole_number(h))];
%     end
end
p = 1;
Pool = {};
tic


for i = 1:length(u_list)
    y = eval(M12E_unit_list.data{u_list(i),4});
    if y.stimulus_ch1(1,3) == dB && y.stimulus_ch1(1,8) == 4000
        
        unit_file_name = 'M12Eu000';
        unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
        x = load(unit_file_name);
        if length(x.s_unit.spiketimes)>10
            if ~strcmp(y.analysis_type,'User: Trill_list.txt')  %temp solution
                
                if isfield(x.s_unit, 'templates') % Check for existing templates
                    if ~isempty(x.s_unit.templates)
                        Pool(p).best_ch = x.s_unit.templates.best_ch;
                    end
                else
                    data = zeros(64,60);
                    for ch = 1:64
                        temp = [];
                        if size(x.s_unit.waveforms,1) == 64
                            temp(:,:) = x.s_unit.waveforms(ch,:,:);
                            
                        else
                            temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
                            
                        end
                        [Y,Sig,X] = svd(temp,'econ');
                        %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                        k = 1:3;
                        P = Y(:,k)*Sig(k,k)*X(:,k)';
                        data(ch,:) = mean(P,2).';
                    end
                    data = abs(data);
                    [~,Pool(p).best_ch] = max(mean(data(:,10:40),2));
                end
                
                %extract waveforms
                if size(x.s_unit.waveforms,1) == 64
                    Pool(p).waveforms(:,:) = x.s_unit.waveforms(Pool(p).best_ch,:,:);
                    
                else
                    Pool(p).waveforms(:,:) = x.s_unit.waveforms{1}(Pool(p).best_ch,:,:);
                end
                %extract spiketimes and other information
                Pool(p).spiketimes = x.s_unit.spiketimes;
                Pool(p).xb = y;
                prev_unit_fn = [unit_file_name(1:end-4) '_prev.mat'];
                if isfield(x.s_unit,'start_times')
                    Pool(p).stim_times = [x.s_unit.start_times x.s_unit.end_times];
                else
                    z = load(prev_unit_fn);
                    Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
                end
                Pool(p).unit_nb = u_list(i);
                
                for ll = 1:length(M12E_neurons_list)
                    if sum(M12E_neurons_list(ll,:) == u_list(i)) ==1
                        Pool(p).neuron_nb = ll;
                    end
                end
                
                
                p = p+1;
                
            end
        end
    end
    toc
end
toc

    
for n = 1:length(Pool)
    if isempty(Pool(n).neuron_nb)
        Pool(n).neuron_nb = -1;
    end
end

%% Determine best frequency



raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_timesSU = {};
rate.stim = {};
rate.pre = {};
rate.PSTH = {};
unique_neuron_list = unique([Pool.neuron_nb]);
neuron_check_list = [];
% figure
average = [];
for p = 1:length(Pool)
    PreStim = Pool(p).xb.pre_stimulus_record_time*1e-3; %s
    PostStim = Pool(p).xb.post_stimulus_record_time*1e-3; %s
    
    stim_info = Pool(p).xb.data(find(Pool(p).xb.data(:,3) == 1 & Pool(p).xb.data(:,4) == -1),:);
    
    
    data_new = stim_info;
    StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
    
    nreps = Pool(p).xb.stimulus_ch1(1,4);
    nStim = max(Pool(p).xb.stimulus_ch1(:,1));
    %             nreps = 10;
    
    %
    % nStim = nStim + max(x2.stimulus_ch1(:,1));
    %
    TotalReps = nStim*nreps;
    false_start = length(Pool(p).stim_times)-TotalReps;
    
    while false_start <0
        nreps = nreps-1;
        TotalReps = nStim*nreps;
        false_start = length(Pool(p).stim_times)-TotalReps;
 
    end
    
    start_stim_times = Pool(p).stim_times(false_start+1:end,1);
    end_stim_times = Pool(p).stim_times(false_start+1:end,2);
    
    
    
    %             for id = 1:length(SU)
    
    
    
    raster.stim{p} = [];
    raster.rep{p} = [];
    raster.spikes{p} = [];
    spikes_pooled{p} = [];
    %                 rate_stim{id} = [];
    %                 if obj.manual == 0
    %                     spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
    %                 else
    %                     spike_timesSU{id} = obj.spike_times{id};
    %                 end
    for rep = 1:TotalReps
        %for raster
        %                     if id ==4
        %                         spike_timesSU{id} = a.spike_times{4};
        %                     end
        spikes1 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep)-PreStim & ...
            Pool(p).spiketimes<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{p} = [spikes_pooled{p} spikes1];
        raster.stim{p} = [raster.stim{p} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{p} = [raster.rep{p} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{p} = [raster.spikes{p} spikes1];
        
        %for rate
        spikes2 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep) & ...
            Pool(p).spiketimes<=end_stim_times(rep))).';
        rate.stim{p}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
        spikes3 = Pool(p).spiketimes(find(Pool(p).spiketimes<=start_stim_times(rep) & ...
            Pool(p).spiketimes>=start_stim_times(rep)-PreStim)).' ;
        rate.pre{p}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
        
        TrialLength = PreStim + StimDur(data_new(rep,1)) + PostStim;
        rep_rate = zeros(1,(round(TrialLength*1e3)));
        spikes4rate = spikes1 + PreStim;
        for st = spikes4rate
            if ceil(st*1e3) <= length(rep_rate)
                rep_rate(1,ceil(st*1e3)) = rep_rate(1,ceil(st*1e3))+1;
            end
            
        end
%         rate_total = [rate_total ; rate*1000];
        % calculate PSTH
        rate.PSTH{p,data_new(rep,1)}(data_new(rep,2),:) = rep_rate*1e3;
        %
    end

    if a_code == 1
        stim_label = Pool(p).xb.stimulus_ch1(:,8);                   
    elseif a_code >2300
        stim_number = Pool(p).xb.stimulus_ch1(:,1);
        stim_length = Pool(p).xb.stimulus_ch1(:,5);
        for i = 9:22
            VT_para(i-8) = length(unique(Pool(p).xb.stimulus_ch1(:,i)));
        end
        stim_label = Pool(p).xb.stimulus_ch1(:,find(VT_para ~= 1)+8);
        if a_code == 2320 || a_code == 2335
            stim_label = Pool(p).xb.stimulus_ch1(:,10);
        end
    elseif a_code >100 && a_code <1000 %user list 
        
        
        
%         if length(unique(Pool(p).xb.stimulus_ch1(:,9)))>length(unique(Pool(p).xb.stimulus_ch1(:,10)))
%             stim_label = Pool(p).xb.stimulus_ch1(:,9);
%         else
%             stim_label = Pool(p).xb.stimulus_ch1(:,10);
%         end
    end

    
    SUrate{p}.mean = mean(rate.stim{p},2) -mean(rate.pre{p},2); %(Spikes/second)
    SUrate{p}.error = std(rate.stim{p},1,2)/sqrt(nreps);
    SUrate{p}.spont = mean(mean(rate.pre{p}));
    
    average = [average SUrate{p}.mean];
%     
    if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
        if ~isempty(Pool(p).neuron_nb)
            errorbar(stim_label,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
            hold on
        end
    end
    plot(stim_label,ones(1,length(stim_label))*SUrate{p}.spont,'--k');
    neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
    drawnow
%     
end

MUA.mean = mean(average,2);
MUA.error = std(average,[],2)/sqrt(length(Pool));

figure
errorbar(stim_label,MUA.mean,MUA.error,'LineWidth',2);
%mean/ multiunit best frequency
max_ind = find(MUA.mean == max(MUA.mean));

BF = stim_label(max_ind);


%% CSD analysis using tjp
addpath('C:\Users\John.Lee\Documents\GitHub\Kilosort-Wanglab\tjp');
clear a;
a = tjp;
pause(0.1);
a.params.session_name = '2019-07-25_13-50-02';
a.params.fpath = fullfile(a.params.fpath1, filesep, a.params.session_name);
a.params.xbz_file_name = 'M12E0427';
% a.initialize;
a.load_data_OE;
pause(0.1);
a.load_LFP;
pause(0.1);
a.CSD_analysis(max_ind);

%% 

unique_neuron_list = unique([Pool.neuron_nb]);
neuron_check_list = [];

for p =1:length(Pool)
    
    if mean(abs(SUrate{p}.mean))> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
        figure
        if ~isempty(Pool(p).neuron_nb)
            errorbar(stim_label,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
            hold on
        end
        plot(stim_label,ones(1,length(stim_label))*SUrate{p}.spont,'--k');
        title(['channel ' num2str(Pool(1,p).best_ch)]);
        
    end
    neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
    drawnow
end

neuron_check_list = [];
% group by channel closeness
%calculate latency at best frequency?
latency = zeros(length(Pool),1);
% i = 0;
for p = 1:length(Pool)
    %         if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
    
    if mean(abs(SUrate{p}.mean))> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
%         i = i+1;
        
        BF = find(abs(SUrate{p}.mean) == max(abs(SUrate{p}.mean)));
        if length(BF)>1
            BF = BF(1);
        end
        spont_raw = rate.PSTH{p,BF}(1:PreStim*1e3);
        spont_mean = mean2(spont_raw);
        spont_std = std2(spont_raw);
        lat_temp = find(mean(abs(rate.PSTH{p,BF}(:,PreStim*1e3:PreStim*1e3+StimDur(BF)*1e3)),1)>spont_mean+spont_std*2);
        if isempty(lat_temp)
            latency(p) = nan;
        else
        latency(p) = lat_temp(1);
        end
    end
        neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];

end
%     latency(p) =





%% Figures 


% distribution of units along the probe%
figure
subplot(1,2,1)
h = histogram(round([(Pool(1,:).best_ch)]/4),[1:16], 'Orientation', 'horizontal');
set(gca, 'YDir','reverse')
subplot(1,2,2)
M = max(max(abs(a.CC{max_ind}))); % abosolute maximum CSD, for the colormap scale
clims = [-M M]; % gives the upper and lower limit for the colormap
% imagesc(a.CC{max_ind});
imagesc(a.CC{max_ind}(:,6000:15000),clims); % CSD as heatmap
colormap(flipud(jet)); % blue = source; red = sink
colorbar('EastOutside');

% average PSTH of units in different areas.
% for p = 1:length(Pool)
%     figure
%     smooth_PSTH{p} = [];
%      xs = 1:600;
%     h = 10;
%         BF_ind = find(SUrate{p}.mean == max(SUrate{p}.mean));
% 
%     for i = 1:600
%         ys(i) = gaussian_kern_reg(xs(i),xs,mean(rate.PSTH{p,BF_ind}),h);
%     end
%     smooth_PSTH{p} = ys;
%     plot(smooth_PSTH{p})
%         title(['channel ' num2str(Pool(1,p).best_ch)]);
% drawnow 
% end


%Onset latency and probe position
figure
scatter(latency,round([Pool(1,:).best_ch]/4));
set(gca, 'YDir','reverse')























        