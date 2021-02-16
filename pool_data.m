
function [Pool, rate, raster] = pool_data(a_code,hole_number,varargin)

%% Comments 
% If a_code is 1, vargargin{1} should be dB attenuation (20,40,60)

animal_name = 'M60F';
addpath(fullfile('D:\DATA', filesep, animal_name, filesep,'Units')); % path to Units
addpath(fullfile('D:\DATA', filesep, animal_name, filesep,'Experiments')); %path to xbz files
load([animal_name '_neurons_list.mat']);

% addpath('D:\DATA\M12E\Units'); % path to units
% addpath('D:\DATA\M12E\Experiments'); %path to xbz files
% addpath('Y:\Units_Backup'); %path to units server
% load('M12E_neurons_list.mat'); 
% note, change neurons_list file from M12E_neurons_list to neurons_list;
% change/ check for M12E_unit_list file file name variable names etc


tic
if ~exist([animal_name '_unit_list_new.mat'])
    load([animal_name '_unit_list.mat']);
    fprintf('Time %3.0fs. adding analysis code... \n', toc);
    
    
    for i  = 1:length(unit_list.data)
        
        x = eval(unit_list.data{i,4});
%         load([unit_list.data{i,4} '.mat']);
        unit_list.data{i,6} = x.analysis_code;
        if strcmp(x.analysis_type,'User: Trill_list.txt')
            unit_list.data{i,6} = 286;
        end
        unit_list.data{i,7} = x.analysis_type;
        unit_list.data{i,8} = x.hole_number;
        unit_list.data{i,9} = x.track_number;
        if mod(i,100) == 0
            fprintf(['Time %3.0fs. adding analysis code...' num2str(i) '/' num2str(length(unit_list.data)) '\n'], toc);
        end
            
    end
    fprintf('Time %3.0fs. adding analysis code... Done! \n', toc);
    
    % fix hole/track numbers
    
    
    
else
    load([animal_name '_unit_list_new.mat']);
end

toc
% for i = 2278: 2472
%    unit_list.data{i,8} = 4;
%     unit_list.data{i,9} = 6;
% end
% 





%% The following code can be used universally, to extract spike and stim times from unit data.
% dB = 40;
% a_code = 1;
% hole_number = 10;

tic
if size(hole_number) == 1
    u_list = find([unit_list.data{:,6}] == a_code & [unit_list.data{:,8}] == hole_number);
else
    u_list = [];
    for h = 1:length(hole_number)
        u_list = [u_list find([unit_list.data{:,6}] == a_code & [unit_list.data{:,8}] == hole_number(h))];
    end
end
p = 1;
Pool = {};
tic


if a_code == 1 || a_code == 62
    dB = varargin{1};
    lens = varargin{2}; %length of stim in ms
    for i = 1:length(u_list)
        y = eval(unit_list.data{u_list(i),4});
        if y.stimulus_ch1(1,5) ==lens && y.stimulus_ch1(1,3) == dB
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
        toc
    end
    toc
    
else

    for i = 1:length(u_list)
        y = eval(unit_list.data{u_list(i),4});
%         if y.stimulus_ch1(1,3) == dB
            unit_file_name = [animal_name 'u00000'];
            unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
            x = load(unit_file_name);
            if length(x.s_unit.spiketimes)>10
%                 if ~strcmp(y.analysis_type,'User: Trill_list.txt')  %temp solution
                    
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
                    prev_unit_fn1 = [unit_file_name(1:end-4) '_prev.mat'];
                    if isfield(x.s_unit,'start_times')
                        Pool(p).stim_times = [x.s_unit.start_times x.s_unit.end_times];
                    else
                        if exist(prev_unit_fn1)
                            z = load(prev_unit_fn1);
                        else
                            
                            prev_unit_fn = [prev_unit_fn1(1:4) prev_unit_fn1(7:end)];
                            if exist(prev_unit_fn)
                                z = load(prev_unit_fn);
                            else
                                prev_unit_fn = [prev_unit_fn1(1:5) prev_unit_fn1(8:end)];
                                z = load(prev_unit_fn);
                                
                            end
                        end
                        Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
                    end
                    Pool(p).unit_nb = u_list(i);
                    
                    Pool(p).neuron_nb = [];
                    for ll = 1:length(neurons_list)
                        if sum(neurons_list(ll,2:end) == u_list(i)) ==1
                            Pool(p).neuron_nb = ll;                            
                        end
                    end
                    Pool(p).hole_nb = unit_list.data{u_list(i),8};
                    Pool(p).track_nb = unit_list.data{u_list(i),9};
                    
                    p = p+1;
                    
                    %                 end
                    %             end
        end
        toc
    end
    toc
end
    
    
    
for n = 1:length(Pool)
    if isempty(Pool(n).neuron_nb)
        Pool(n).neuron_nb = -1;
    end
end
        
                

%% Raster and Tuning curves

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
    
%     if a_code == 1
%         stim_label = Pool(p).xb.stimulus_ch1(:,8);                   
%     elseif a_code >2300
%         stim_number = Pool(p).xb.stimulus_ch1(:,1);
%         stim_length = Pool(p).xb.stimulus_ch1(:,5);
%         for i = 9:22
%             VT_para(i-8) = length(unique(Pool(p).xb.stimulus_ch1(:,i)));
%         end
%         stim_label = Pool(p).xb.stimulus_ch1(:,find(VT_para ~= 1)+8);
%         if a_code == 2320 || a_code == 2335
%             stim_label = Pool(p).xb.stimulus_ch1(:,10);
%         end
%     elseif a_code >100 && a_code <1000 %user list 
%         
%         
%         
% %         if length(unique(Pool(p).xb.stimulus_ch1(:,9)))>length(unique(Pool(p).xb.stimulus_ch1(:,10)))
% %             stim_label = Pool(p).xb.stimulus_ch1(:,9);
% %         else
% %             stim_label = Pool(p).xb.stimulus_ch1(:,10);
% %         end
%     end

    
%     SUrate{p}.mean = mean(rate.stim{p},2) -mean(rate.pre{p},2); %(Spikes/second)
%     SUrate{p}.error = std(rate.stim{p},1,2)/sqrt(nreps);
%     SUrate{p}.spont = mean(mean(rate.pre{p}));
    
    
%     if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
%         if ~isempty(Pool(p).neuron_nb)
%             errorbar(stim_label,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
%             hold on
%         end
%     end
    % plot(stim_label,ones(1,length(stim_label))*SUrate{p}.spont,'--k');
    neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
%     drawnow
    
end












% 
% 
% 
% 
% 
% 
% 
% 
% 
% tic
% if ~exist('M12E_unit_list_new.mat')
%     load('M12E_unit_list.mat');
%     fprintf('Time %3.0fs. adding analysis code... \n', toc);
%     
%     
%     for i  = 1:length(M12E_unit_list.data)
%         
%         x = eval(M12E_unit_list.data{i,4});
%         load([ M12E_unit_list.data{i,4} '.mat']);
%         M12E_unit_list.data{i,6} = x.analysis_code;
%         if strcmp(x.analysis_type,'User: Trill_list.txt')
%             M12E_unit_list.data{i,6} = 286;
%         end
%         M12E_unit_list.data{i,7} = x.analysis_type;
%         M12E_unit_list.data{i,8} = x.hole_number;
%         M12E_unit_list.data{i,9} = x.track_number;
%         if mod(i,100) == 0
%             fprintf(['Time %3.0fs. adding analysis code...' num2str(i) '/' num2str(length(M12E_unit_list.data)) '\n'], toc);
%         end
%             
%     end
%     fprintf('Time %3.0fs. adding analysis code... Done! \n', toc);
%     
% else
%     load('M12E_unit_list_new.mat');
% end
% 
% 
% 
% 
% 
% 
% 
% %% The following code can be used universally, to extract spike and stim times from unit data.
% % dB = 40;
% % a_code = 1;
% % hole_number = 10;
% 
% tic
% if size(hole_number) == 1
%     u_list = find([M12E_unit_list.data{:,6}] == a_code & [M12E_unit_list.data{:,8}] == hole_number);
% else
%     u_list = [];
%     for h = 1:length(hole_number)
%         u_list = [u_list find([M12E_unit_list.data{:,6}] == a_code & [M12E_unit_list.data{:,8}] == hole_number(h))];
%     end
% end
% p = 1;
% Pool = {};
% tic
% 
% 
% if a_code == 1 || a_code == 62
%     dB = varargin{1};
%     
%     for i = 1:length(u_list)
%         y = eval(M12E_unit_list.data{u_list(i),4});
%         if y.stimulus_ch1(1,3) == dB
%             unit_file_name = 'M12Eu00000';
%             unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
%             x = load(unit_file_name);
%             if length(x.s_unit.spiketimes)>10
%                 if ~strcmp(y.analysis_type,'User: Trill_list.txt')  %temp solution
%                     
%                     if isfield(x.s_unit, 'templates') % Check for existing templates
%                         if ~isempty(x.s_unit.templates)
%                             Pool(p).best_ch = x.s_unit.templates.best_ch;
%                         end
%                     else
%                         data = zeros(64,60);
%                         for ch = 1:64
%                             temp = [];
%                             if size(x.s_unit.waveforms,1) == 64
%                                 temp(:,:) = x.s_unit.waveforms(ch,:,:);
%                                 
%                             else
%                                 temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
%                                 
%                             end
%                             [Y,Sig,X] = svd(temp,'econ');
%                             %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
%                             k = 1:3;
%                             P = Y(:,k)*Sig(k,k)*X(:,k)';
%                             data(ch,:) = mean(P,2).';
%                         end
%                         data = abs(data);
%                         [~,Pool(p).best_ch] = max(mean(data(:,10:40),2));
%                     end
%                     
%                     %extract waveforms
%                     if size(x.s_unit.waveforms,1) == 64
%                         Pool(p).waveforms(:,:) = x.s_unit.waveforms(Pool(p).best_ch,:,:);
%                         
%                     else
%                         Pool(p).waveforms(:,:) = x.s_unit.waveforms{1}(Pool(p).best_ch,:,:);
%                     end
%                     %extract spiketimes and other information
%                     Pool(p).spiketimes = x.s_unit.spiketimes;
%                     Pool(p).xb = y;
%                     
%                     prev_unit_fn = [unit_file_name(1:end-4) '_prev.mat'];
%                     if isfield(x.s_unit,'start_times')
%                         Pool(p).stim_times = [x.s_unit.start_times x.s_unit.end_times];
%                     else
%                         if exist(prev_unit_fn)
%                             z = load(prev_unit_fn);
%                         else
%                             if strcmp(prev_unit_fn(7),'0')
%                                 prev_unit_fn = [prev_unit_fn(1:5) prev_unit_fn(8:end)];
%                             else
%                                 prev_unit_fn = [prev_unit_fn(1:4) prev_unit_fn(7:end)];
%                             end
%                             z = load(prev_unit_fn);
%                             
%                         end
%                         z = load(prev_unit_fn);
%                         Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
%                     end
%                     Pool(p).unit_nb = u_list(i);
%                     Pool(p).neuron_nb = [];
%                     for ll = 1:length(M12E_neurons_list)
%                         if sum(M12E_neurons_list(ll,2:end) == u_list(i)) ==1
%                             Pool(p).neuron_nb = ll;                            
%                         end
%                     end
% %                     
%                     Pool(p).hole_nb = M12E_unit_list.data{u_list(i),8};
%                     Pool(p).track_nb = M12E_unit_list.data{u_list(i),9};
%                     
%                     p = p+1;
%                     
%                 end
%             end
%         end
%         toc
%     end
%     toc
%     
% else
% 
%     for i = 1:length(u_list)
%         y = eval(M12E_unit_list.data{u_list(i),4});
% %         if y.stimulus_ch1(1,3) == dB
%             unit_file_name = 'M12Eu00000';
%             unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
%             x = load(unit_file_name);
%             if length(x.s_unit.spiketimes)>10
% %                 if ~strcmp(y.analysis_type,'User: Trill_list.txt')  %temp solution
%                     
%                     if isfield(x.s_unit, 'templates') % Check for existing templates
%                         if ~isempty(x.s_unit.templates)
%                             Pool(p).best_ch = x.s_unit.templates.best_ch;
%                         end
%                     else
%                         data = zeros(64,60);
%                         for ch = 1:64
%                             temp = [];
%                             if size(x.s_unit.waveforms,1) == 64
%                                 temp(:,:) = x.s_unit.waveforms(ch,:,:);
%                                 
%                             else
%                                 temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
%                                 
%                             end
%                             [Y,Sig,X] = svd(temp,'econ');
%                             %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
%                             k = 1:3;
%                             P = Y(:,k)*Sig(k,k)*X(:,k)';
%                             data(ch,:) = mean(P,2).';
%                         end
%                         data = abs(data);
%                         [~,Pool(p).best_ch] = max(mean(data(:,10:40),2));
%                     end
%                     
%                     %extract waveforms
%                     if size(x.s_unit.waveforms,1) == 64
%                         Pool(p).waveforms(:,:) = x.s_unit.waveforms(Pool(p).best_ch,:,:);
%                         
%                     else
%                         Pool(p).waveforms(:,:) = x.s_unit.waveforms{1}(Pool(p).best_ch,:,:);
%                     end
%                     %extract spiketimes and other information
%                     Pool(p).spiketimes = x.s_unit.spiketimes;
%                     Pool(p).xb = y;
%                     prev_unit_fn1 = [unit_file_name(1:end-4) '_prev.mat'];
%                     if isfield(x.s_unit,'start_times')
%                         Pool(p).stim_times = [x.s_unit.start_times x.s_unit.end_times];
%                     else
%                         if exist(prev_unit_fn1)
%                             z = load(prev_unit_fn1);
%                         else
%                             
%                             prev_unit_fn = [prev_unit_fn1(1:4) prev_unit_fn1(7:end)];
%                             if exist(prev_unit_fn)
%                                 z = load(prev_unit_fn);
%                             else
%                                 prev_unit_fn = [prev_unit_fn1(1:5) prev_unit_fn1(8:end)];
%                                 z = load(prev_unit_fn);
%                                 
%                             end
%                         end
%                         Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
%                     end
%                     Pool(p).unit_nb = u_list(i);
%                     
%                     Pool(p).neuron_nb = [];
%                     for ll = 1:length(M12E_neurons_list)
%                         if sum(M12E_neurons_list(ll,2:end) == u_list(i)) ==1
%                             Pool(p).neuron_nb = ll;                            
%                         end
%                     end
%                     Pool(p).hole_nb = M12E_unit_list.data{u_list(i),8};
%                     Pool(p).track_nb = M12E_unit_list.data{u_list(i),9};
%                     
%                     p = p+1;
%                     
%                     %                 end
%                     %             end
%         end
%         toc
%     end
%     toc
% end
%     
%     
%     
% for n = 1:length(Pool)
%     if isempty(Pool(n).neuron_nb)
%         Pool(n).neuron_nb = -1;
%     end
% end
%         
%                 
% 
% %% Raster and Tuning curves
% 
% % p = 10;
% rate = {};
% raster = {};
% raster.stim = {};
% raster.rep = {};
% raster.spikes = {};
% spikes_pooled = {};
% spike_timesSU = {};
% rate.stim = {};
% rate.pre = {};
% rate.post = {};
% rate.PSTH = {};
% 
% unique_neuron_list = unique([Pool.neuron_nb]);
% neuron_check_list = [];
% % figure
% for p = 1:length(Pool)
%     PreStim = Pool(p).xb.pre_stimulus_record_time*1e-3; %s
%     PostStim = Pool(p).xb.post_stimulus_record_time*1e-3; %s
%     
%     stim_info = Pool(p).xb.data(find(Pool(p).xb.data(:,3) == 1 & Pool(p).xb.data(:,4) == -1),:);
%     
%     
%     data_new = stim_info;
%     StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
%     
%     nreps = Pool(p).xb.stimulus_ch1(1,4);
%     nStim = max(Pool(p).xb.stimulus_ch1(:,1));
%     %             nreps = 10;
%     
%     %
%     % nStim = nStim + max(x2.stimulus_ch1(:,1));
%     %
%     TotalReps = nStim*nreps;
%     false_start = length(Pool(p).stim_times)-TotalReps;
%     
%     while false_start <0
%         nreps = nreps-1;
%         TotalReps = nStim*nreps;
%         false_start = length(Pool(p).stim_times)-TotalReps;
%  
%     end
%     
%     start_stim_times = Pool(p).stim_times(false_start+1:end,1);
%     end_stim_times = Pool(p).stim_times(false_start+1:end,2);
%     
%     
%     
%     %             for id = 1:length(SU)
%     
%     
%     
%     raster.stim{p} = [];
%     raster.rep{p} = [];
%     raster.spikes{p} = [];
%     spikes_pooled{p} = [];
%     %                 rate_stim{id} = [];
%     %                 if obj.manual == 0
%     %                     spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
%     %                 else
%     %                     spike_timesSU{id} = obj.spike_times{id};
%     %                 end
%     for rep = 1:TotalReps
%         %for raster
%         %                     if id ==4
%         %                         spike_timesSU{id} = a.spike_times{4};
%         %                     end
%         spikes1 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep)-PreStim & ...
%             Pool(p).spiketimes<=end_stim_times(rep)+ PostStim)).';
%         spikes1 = spikes1 - start_stim_times(rep);
%         spikes_pooled{p} = [spikes_pooled{p} spikes1];
%         raster.stim{p} = [raster.stim{p} data_new(rep,1)*ones(size(spikes1))];
%         raster.rep{p} = [raster.rep{p} data_new(rep,2)*ones(size(spikes1))];
%         raster.spikes{p} = [raster.spikes{p} spikes1];
%         
%         %for rate
%         spikes2 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep) & ...
%             Pool(p).spiketimes<=end_stim_times(rep))).';
%         rate.stim{p}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
%         spikes3 = Pool(p).spiketimes(find(Pool(p).spiketimes<=start_stim_times(rep) & ...
%             Pool(p).spiketimes>=start_stim_times(rep)-PreStim)).' ;
%         spikes4 = Pool(p).spiketimes(find(Pool(p).spiketimes<=end_stim_times(rep)+ PostStim & ...
%             Pool(p).spiketimes>=end_stim_times(rep))).' ;
%         
%         rate.pre{p}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
%         rate.post{p}(data_new(rep,1),data_new(rep,2)) = length(spikes4)/PostStim;
%         
%         TrialLength = PreStim + StimDur(data_new(rep,1)) + PostStim;
%         rep_rate = zeros(1,(round(TrialLength*1e3)));
%         spikes4rate = spikes1 + PreStim;
%         for st = spikes4rate
%             if st<0
%                 st1 = 1;
%                 if ceil(st1*1e3) <= length(rep_rate)
%                     rep_rate(1,ceil(st1*1e3)) = rep_rate(1,ceil(st1*1e3))+1;
%                 end
%             
%             elseif ceil(st*1e3) <= length(rep_rate)
%                 rep_rate(1,ceil(st*1e3)) = rep_rate(1,ceil(st*1e3))+1;
%             end
%             
%         end
% %         rate_total = [rate_total ; rate*1000];
%         % calculate PSTH
%         rate.PSTH{p,data_new(rep,1)}(data_new(rep,2),:) = rep_rate*1e3;
%         %
%         
%         
%     end
%     
% %     if a_code == 1
% %         stim_label = Pool(p).xb.stimulus_ch1(:,8);                   
% %     elseif a_code >2300
% %         stim_number = Pool(p).xb.stimulus_ch1(:,1);
% %         stim_length = Pool(p).xb.stimulus_ch1(:,5);
% %         for i = 9:22
% %             VT_para(i-8) = length(unique(Pool(p).xb.stimulus_ch1(:,i)));
% %         end
% %         stim_label = Pool(p).xb.stimulus_ch1(:,find(VT_para ~= 1)+8);
% %         if a_code == 2320 || a_code == 2335
% %             stim_label = Pool(p).xb.stimulus_ch1(:,10);
% %         end
% %     elseif a_code >100 && a_code <1000 %user list 
% %         
% %         
% %         
% % %         if length(unique(Pool(p).xb.stimulus_ch1(:,9)))>length(unique(Pool(p).xb.stimulus_ch1(:,10)))
% % %             stim_label = Pool(p).xb.stimulus_ch1(:,9);
% % %         else
% % %             stim_label = Pool(p).xb.stimulus_ch1(:,10);
% % %         end
% %     end
% 
%     
% %     SUrate{p}.mean = mean(rate.stim{p},2) -mean(rate.pre{p},2); %(Spikes/second)
% %     SUrate{p}.error = std(rate.stim{p},1,2)/sqrt(nreps);
% %     SUrate{p}.spont = mean(mean(rate.pre{p}));
%     
%     
% %     if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
% %         if ~isempty(Pool(p).neuron_nb)
% %             errorbar(stim_label,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
% %             hold on
% %         end
% %     end
%     % plot(stim_label,ones(1,length(stim_label))*SUrate{p}.spont,'--k');
%     neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
% %     drawnow
%     
% end
% 


