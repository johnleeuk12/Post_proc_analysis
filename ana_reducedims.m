function [projs mse like lat] = ana_reducedims(D, dims)  

% D = DD;
% dims = 1:10;

% remove low firing rate neurons
mean_thresh= 1;
m = mean([D.data],2)*1e3;
keep_neurons = m >= mean_thresh;


% Remove low firing rate neurons

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



if (binWidth ~= 1 || use_sqrt)
    [d(1:length(D)).spikes] = deal(D.data);
    for itrial = 1:length(D)
        d(itrial).trialId = itrial;
    end
    
    s = getSeq(d, binWidth, 'useSqrt', use_sqrt);
    
    [D.data] = deal(s.y);
end



% Trial-average neural trajectories
D = trial_average(D);





%
% Smooth data if necessary (automatically zero if GPFA selected)
%     if get(handles.kern_slider,'Value') ~= 0

smooth_bin = 10; % ['Smoothing kernel width: '  num2str(value) 'ms']);

if smooth_bin ~= 1
    for i = 1:length(D)
        D(i).data = smoother(D(i).data, smooth_bin , binWidth);
    end
end
%
lat   = [];

like_fold = 0;
mse_fold = 0;  % these store the likelihood and mse for one fold
like  = zeros(1,length(dims));
mse   = zeros(1,length(dims)); % keeps a running total of likelihood and mse
projs = cell(1,length(dims));


% prepare cross-validation technique
% break up into folds
cv_trials = randperm(length(D));
mask = false(1, length(D));
fold_indices = floor(linspace(1,length(D)+1, 4));  %splits it up into three chunks

for idim = 1 :length(dims)
    for ifold = 1:3 % three-fold cross-validation
        fprintf(['Cross validating... dim ' ...
            num2str(dims(idim)) ' fold ' num2str(ifold) '\n'])
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask = mask;
        test_mask(cv_trials(fold_indices(ifold):fold_indices(ifold+1)-1)) = true;
        train_mask = ~test_mask;
        [mse_fold p] = PCACV(D, dims(idim), ifold, train_mask, test_mask);
        
        
        % add up the likelihood and LNO errors across folds
        mse(idim) = mse(idim) + mse_fold;
        like(idim) = like(idim) + like_fold;
        
        
        
        if (ifold == 1)
            projs{idim} = p;
        end
    end
end

[u sc lat] = pca([D.data]');
end


% -------- Helper Functions
function [mse, projs, lat] = PCACV(D, dim, fold, train_mask, test_mask)

    [train_data, test_data, forProj] = prepare_cv_data(D, train_mask, test_mask);
    [u, sc, lat] = pca(train_data');
    params.L = u(:,1:dim);
    params.d = mean(train_data,2);
    projs = [];
    
    if (fold == 1 && dim > 1) % keep the projections of first fold
        % All of the data projected into low-D space
        allprojs = sc(:,1:dim)';

        if (isstruct(forProj)) % trajectories
            index = 1;
            for itrial=1:length(forProj)
                projs(itrial).data = allprojs(:,index:index+size(forProj(itrial).data,2)-1);
                index = index + size(forProj(itrial).data,2);
            end
            [projs(1:end).type] = deal('traj');
        else   %clusters
            projs(1).data = allprojs;
            projs(1).type = 'state';
        end
    end

    cvdata = cosmoother_pca(test_data,params);
    mse = sum(sum((cvdata-test_data).^2));

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

function [train_data, test_data, forProj] = prepare_cv_data(D, train_mask, test_mask)
% helper function to prepare the cv data
% finds the train and test data, as well as the struct for forProjs
% This is useful to handle both trajs and states

    train_data = [];
    test_data = [];
%     if (strcmp(D(1).type, 'traj')) % data is trajectories
        train_data = [D(train_mask).data];
        test_data = [D(test_mask).data];
        forProj = D(train_mask);
%     else  % data is clusters
%         data = [D.data];
%         train_data = data(:,train_mask);
%         test_data = data(:,test_mask);
%         forProj = data(:,train_mask);
%     end
end
            
            