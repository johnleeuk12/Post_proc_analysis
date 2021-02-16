function [Pool, rate, raster] =  Pool_BF(hole_number)
%%
%{
2/15/2021
determining BF and amplitude in response to PT and BP noise

%}
%% Loading data
% run Pool_data for unit_list_new first
animal_name = 'M60F';
addpath(fullfile('D:\DATA', filesep, animal_name, filesep,'Units')); % path to Units
addpath(fullfile('D:\DATA', filesep, animal_name, filesep,'Experiments')); %path to xbz files
load([animal_name '_neurons_list.mat']);
load([animal_name '_unit_list_new.mat']);

% for i = 2278: 2472
%    unit_list.data{i,8} = 4;
%     unit_list.data{i,9} = 6;
% end
%

%%

% hole_number = 4;
a_code = [1 62];
u_list = [];

tic
fprintf('determining u_list...')
if size(hole_number) == 1
    u_list = find([unit_list.data{:,6}] == a_code(1) & [unit_list.data{:,8}] == hole_number);
    u_list = [u_list find([unit_list.data{:,6}] == a_code(2) & [unit_list.data{:,8}] == hole_number)];
else
    u_list = [];
    for h = 1:length(hole_number)
        u_list = [u_list find([unit_list.data{:,6}] == a_code & [unit_list.data{:,8}] == hole_number(h))];
    end
end

fprintf('done!. elapsed time is: %.2f seconds. \n',toc')


p = 1;
Pool = {};

for i = 1:length(u_list)
    if mod(i,100) ==1
        fprintf(['%4d /' num2str(length(u_list)) ' time : %6.2f sec \n'],i,toc')
    end
    y = eval(unit_list.data{u_list(i),4});
    %     if y.stimulus_ch1(1,5) ==lens && y.stimulus_ch1(1,3) == dB
    unit_file_name = [animal_name 'u00000'];
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
                if exist(prev_unit_fn)
                    z = load(prev_unit_fn);
                else
                    if strcmp(prev_unit_fn(7),'0')
                        prev_unit_fn = [prev_unit_fn(1:5) prev_unit_fn(8:end)];
                    else
                        prev_unit_fn = [prev_unit_fn(1:4) prev_unit_fn(7:end)];
                    end
                    z = load(prev_unit_fn);
                    
                end
                z = load(prev_unit_fn);
                Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
            end
            Pool(p).unit_nb = u_list(i);
            Pool(p).neuron_nb = [];
            for ll = 1:length(neurons_list)
                if sum(neurons_list(ll,2:end) == u_list(i)) ==1
                    Pool(p).neuron_nb = ll;
                end
            end
            %
            Pool(p).hole_nb = unit_list.data{u_list(i),8};
            Pool(p).track_nb = unit_list.data{u_list(i),9};
            
            p = p+1;
            
        end
    end
    
    
end

%% Raster and Tuning curves

fprintf('Rate and Raster. time : %4.2f sec. \n',toc')



% p = 10;
rate = {};
raster = {};
raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_timesSU = {};
rate.stim = {};
rate.pre = {};
rate.post = {};
rate.PSTH = {};

unique_neuron_list = unique([Pool.neuron_nb]);
neuron_check_list = [];
% figure
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
        spikes4 = Pool(p).spiketimes(find(Pool(p).spiketimes<=end_stim_times(rep)+ PostStim & ...
            Pool(p).spiketimes>=end_stim_times(rep))).' ;
        
        rate.pre{p}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
        rate.post{p}(data_new(rep,1),data_new(rep,2)) = length(spikes4)/PostStim;
        
        TrialLength = PreStim + StimDur(data_new(rep,1)) + PostStim;
        rep_rate = zeros(1,(round(TrialLength*1e3)));
        spikes4rate = spikes1 + PreStim;
        for st = spikes4rate
            if st<0
                st1 = 1;
                if ceil(st1*1e3) <= length(rep_rate)
                    rep_rate(1,ceil(st1*1e3)) = rep_rate(1,ceil(st1*1e3))+1;
                end
            
            elseif ceil(st*1e3) <= length(rep_rate)
                rep_rate(1,ceil(st*1e3)) = rep_rate(1,ceil(st*1e3))+1;
            end
            
        end
%         rate_total = [rate_total ; rate*1000];
        % calculate PSTH
        rate.PSTH{p,data_new(rep,1)}(data_new(rep,2),:) = rep_rate*1e3;
        %
        
        
    end

    neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
%     drawnow
    
end

for n = 1:length(Pool)
    if isempty(Pool(n).neuron_nb)
        Pool(n).neuron_nb = -1;
    end
end
     
fprintf('done! time: %4.2f sec. \n',toc')






