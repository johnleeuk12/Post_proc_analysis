function pool_C_BF()


%% calculate BF for best_ind or for 40db, per unit


N_list = out.data{1,1}(:,1);

animal_name = 'M56E';

% savedir = fullfile('D:\DATA',filesep,animal_name,filesep,'Analysis\ana_tones\data');
% addpath(savedir);
% load([savedir '\neurons_loc_tag.mat']);
% for l = 1:length(neurons_loc_tag)
%     if isempty(neurons_loc_tag(l).neuron_nb)
%         neurons_loc_tag(l).neuron_nb = -1;
%         neurons_loc_tag(l).hole_nb = -1;
%         neurons_loc_tag(l).track_nb = -1;
%     end
% end


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


C = {};
C.pool.neurons_info = out.data{1,1};
C.pool.data = B;
C.pool.data_tags = {'db_ind'...
                    'BF, kHz'...
                    'peak_latency, ms'...
                    'peak_FR, spks/s'...
                    'spont_FR, spks/s'...
                    'min_latency, ms'...
                    'nb of db'};
