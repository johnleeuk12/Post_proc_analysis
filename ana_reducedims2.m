function [newD, C, lat,explained,keep_neurons] = ana_reducedims2(D, dims,smooth_bin,low_f)  
tic

addpath('C:\Users\John.Lee\Documents\GitHub\DataHigh');
addpath('C:\Users\John.Lee\Documents\GitHub\DataHigh\gpfa');
addpath('util');
% D = DD;
% dims = 1:10;

% remove low firing rate neurons

fprintf(['removing low firing rate neurons ... time : %6.2f sec \n'],toc')

mean_thresh= low_f;
m = mean([D.data],2)*1e3;
keep_neurons = m >= mean_thresh;



for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(keep_neurons,:);
end


% find minimum trial length
min_trial_length = inf;
for itrial = 1:length(D)
    if (size(D(itrial).data,2) < min_trial_length)
        min_trial_length = size(D(itrial).data,2);
    end
end

binWidth =min(20, min_trial_length);
use_sqrt = 0;

fprintf(['gathering spikes ... time : %6.2f sec \n'],toc')

% gather spikes in 20ms bins. 
if (binWidth ~= 1 || use_sqrt)
    [d(1:length(D)).spikes] = deal(D.data);
    for itrial = 1:length(D)
        d(itrial).trialId = itrial;
    end
    
    s = getSeq(d, binWidth, 'useSqrt', use_sqrt);
    
    [D.data] = deal(s.y);
end


fprintf(['trial average and normalizing ... time : %6.2f sec \n'],toc')

% Trial-average neural trajectories
D = trial_average(D);

D = normalize_D(D);
% Smooth data if necessary (automatically zero if GPFA selected)
%     if get(handles.kern_slider,'Value') ~= 0

% smooth_bin = 20; % ['Smoothing kernel width: '  num2str(value) 'ms']);

if smooth_bin ~= 1
    for i = 1:length(D)
        D(i).data = smoother(D(i).data, smooth_bin , binWidth);
    end
end
%

% [newD, C, lat,explained] = PCAreduce(D,dims);
[newD, C, lat,explained] = GPFAreduce(D,dims, binWidth);

fprintf(['reduced dimensions ... time : %6.2f sec \n'],toc')


end

%% Helper functions


function [newD, C, lat, explained] = PCAreduce(D,dims)
% PCAREDUCE Internal function for PCA
%   PCAREDUCE(D,DIMS) returns a structure of the same form as D, except
%   the data has been reduced with PCA. All conditions and trials are
%   considered together to get the best joint reduction.

    % Agglomerate all of the conditions, and perform PCA
    alldata = [D.data];
    [u, sc, lat,~,explained] = pca(alldata');

    % For each condition, store the reduced version of each data vector
    index = 0;
    for i=1:length(D)
        D(i).data = sc(index + (1:size(D(i).data,2)),1:dims)';
        index = index + size(D(i).data,2);
    end
    newD = D;
    C = u(:,1:dims);
    lat = cumsum(lat(1:dims)) ./ sum(lat(1:dims));  % eigenvalues
end

function [newD, C, lat, explained] = GPFAreduce(D,dims, binWidth)
% GPFAREDUCE Internal function for GPFA
%   GPFAreduce(D,DIMS) returns a structure of the same form as D, except
%   the data has been reduced with GPFA. All conditions and trials are
%   considered together to get the best joint reduction. Code is based off
%   of GPFAENGINE, by Byron Yu and John Cunningham.

    wb = waitbar(0, 'Performing dim reduction with GPFA');
    [params, newD] = gpfa_engineDH(D,dims,'dimreduce_wb',wb, 'binWidth', binWidth);
    delete(wb);
    C = params.Corth;
    [~, lat, explained] = pcacov(params.C * params.C');
    lat = cumsum(lat(1:dims))./sum(lat(1:dims));
end

function newD = normalize_D(D)

% % soft normalization after 4/5/22
N = size(D(1).data,1);
for n = 1:N
    
    T = zeros(length(D),size(D(1).data,2));
    for c = 1:length(D)
        T(c,:) = D(c).data(n,:);
    end
    
    T = T-max(max(T));
    T = T-mean(T,1);
    for c = 1:length(D)
        D(c).data(n,:) = T(c,:);
    end
    
    
end
newD = D;

%     % before 4/5/2022
%     % removing common responses to enhance difference
%     N = size(D(1).data,1);
% 
%     for n = 1:N
%         T = zeros(length(D),size(D(1).data,2));
%         for c = 1:length(D)
%             T(c,:) = D(c).data(n,:);
%         end
% 
%         for c = 1:length(D)
%             
%             % without soft normalization
% %             D(c).data(n,:) = D(c).data(n,:)-mean(T,1); 
%             % with normalization
%             D(c).data(n,:) = (D(c).data(n,:)-mean(T,1))/max(max(abs(T-mean(T,1)))); 
%         end
%     end
% 
%     newD = D;
end

function newD = trial_average(D)
% Helper function
% Compute trial-averaged neural trajectories

        conditions = unique({D.condition});  % D.condition should exist
                                    % if not given by user, it makes each
                                    % traj its own condition
        if (length(conditions) == length(D)) % each traj is its own cond
              % this means user input trajs without defining cond
              % if every traj had its own condition, the same dim red
              % results would come from single-trials...
              % so make all trials same cond (which is probably what they
              % want to do)
              [D.condition] = deal('1');
              conditions = {'1'};
        end
              
        newD = [];
        for icond = 1:length(conditions)

            % find smallest trial length
            trial_lengths = [];
            trials = find(ismember({D.condition}, conditions{icond}));
            % need to keep track of which trials are in the condition
            
            for itrial = 1:length(trials)
                trial_lengths = [trial_lengths size(D(trials(itrial)).data,2)];
            end
            smallest_trial_length = min(trial_lengths);

            % truncate all trials to smallest trial length
            for itrial = 1:length(trials)
                D(trials(itrial)).data = D(trials(itrial)).data(:,1:smallest_trial_length);
            end

            % average across trials
            m = zeros(size(D(trials(1)).data));
            for itrial = 1:length(trials)
                m = m + D(trials(itrial)).data;
            end
            m = m ./ length(trials);
            
            % place into new struct
            newD(icond).data = m;
            %             newD(icond).type = D(trials(1)).type;
            newD(icond).type = 'traj';
            newD(icond).condition = D(trials(1)).condition;
            newD(icond).epochStarts = D(trials(1)).epochStarts;
            newD(icond).epochColors = D(trials(1)).epochColors;
        end

end

