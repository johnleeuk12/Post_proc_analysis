function Pool_BF()
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

hole_number = 4;
a_code = [1 62];
u_list = [];
if size(hole_number) == 1
    u_list = find([unit_list.data{:,6}] == a_code(1) & [unit_list.data{:,8}] == hole_number);
    u_list = [u_list find([unit_list.data{:,6}] == a_code(2) & [unit_list.data{:,8}] == hole_number)];
else
    u_list = [];
    for h = 1:length(hole_number)
        u_list = [u_list find([unit_list.data{:,6}] == a_code & [unit_list.data{:,8}] == hole_number(h))];
    end
end
p = 1;
Pool = {};

for i = 1:length(u_list)
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
    
    toc
end
toc