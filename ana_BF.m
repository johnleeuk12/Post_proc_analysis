function out = ana_BF(H_nb_list, figure_on,animal_name)
%{

01//04/2022 JHL
output format:
for each dB level, out.data has 5 columns:
the first column is BF, in kHz.
the second column is BF peak latency, in ms.
the third column is BF peak FR, in ms.
the fourth column is BF min latency, in ms.
the fifth column is trial averaged PSTH.
the sixth column is the stimulus set.


%}
















%% CF analysis 03/07/2021
% load('M60F_unit_list_new.mat');
% animal_name = 'M160E';
% H_nb = 4;
% db_ind = 2;4 is 20db % 3 is 40db, 2 is 60db 1 is 80 db

% edges = logspace(3,4.5,25);
% new_edges = zeros(1,39);
% for e = 1:39
%     new_edges(e) = (edges(e) + edges(e+1))/2;
% end



savedir = fullfile('E:\DATA',filesep,animal_name,filesep,'ana_tones\data');
addpath(savedir);
load([savedir '\neurons_loc_tag.mat']);
for l = 1:length(neurons_loc_tag)
    if isempty(neurons_loc_tag(l).neuron_nb)
        neurons_loc_tag(l).neuron_nb = -1;
        neurons_loc_tag(l).hole_nb = -1;
        neurons_loc_tag(l).track_nb = -1;
    end
end

%%
out.tags = {{'nid'},{'80db'},{'60db'},{'40db'},{'20db'},{'Hole and Track nb'}};
out.data = {};
out.data{1} = [];
tuning = {};
N_list = [neurons_loc_tag.neuron_nb];



tic
db_lens = 0;

for H_nb = H_nb_list
    hole_ind = find([neurons_loc_tag.hole_nb] == H_nb);
    tracks = unique([neurons_loc_tag(hole_ind).track_nb]);

    % figure
    for t = 1:length(tracks)
        fprintf(['processing H' num2str(H_nb) 'T' num2str(tracks(t)) '     time : %6.2f sec \n'],toc')
        
        sub_n_list = N_list(find([neurons_loc_tag.track_nb] ==tracks(t) &[neurons_loc_tag.hole_nb] == H_nb ));
        for n =1:length(sub_n_list)
            nid = sub_n_list(n);
            PTname= 'PTn0000';
            PTname = [PTname(1:end-length(num2str(nid))) num2str(nid) '.mat'];
            try
                load(PTname);
                try
                    ana_PT = filename; %currently loaded file is filename. change afterwards
                catch
                end
                db_count = 0; % counter to see if there's at least one non-empty BF_ind
                out_vec = {};
                out_stim = {};
                peak_lat = {}; % peak_latency;
                min_lat = {}; % minimum latency;
                peak_FR = {};
                spont_FR = {};
                % figures
                if figure_on ==  1
                    fi = figure;
                    set(fi, 'Position', [100 200 2400 660]);
                    sgtitle(['unit' num2str(nid) ' H' num2str(H_nb) 'T' num2str(tracks(t))])
                end
                % figures end
                db_count2 = 0;
                %         db_count3 = 0;
                for db_ind1 = 1:4
                    out_vec{db_ind1} = [];
                    tuning{db_ind1} = [];
                    out_stim{db_ind1} = [];
                    peak_lat{db_ind1} = [];
                    min_lat{db_ind1} = [];
                    
                    try
                        if ~isempty(ana_PT{1,db_ind1})
                            try
                                if ~isempty(ana_PT{2,db_ind1})
                                    stim_freq1 = ana_PT{1,db_ind1}.xb.stimulus_ch1(:,8);
                                    stim_freq2 = ana_PT{2,db_ind1}.xb.stimulus_ch1(:,8);
                                    
                                    if length(stim_freq2)> length(stim_freq1)
                                        db_ind2 = 2;
                                    else
                                        db_ind2 = 1;
                                    end
                                else
                                    db_ind2 = 1;
                                end
                            catch
                                db_ind2 = 1;
                            end
                            
                            db_count2 = 1;
                            stim_freq = ana_PT{db_ind2,db_ind1}.xb.stimulus_ch1(:,8);
                            
                            X = [];
                            for st = 1:length(stim_freq)
                                %                         X = [X; mean(ana_PT{db_ind}.PSTH{st},1)-mean(ana_PT{db_ind}.spont)];
                                X = [X; mean(ana_PT{db_ind2,db_ind1}.PSTH{st},1)]; %-mean(ana_PT{db_ind}.spont)]
                                
                            end
                            
                            %                     spont_std = std2(X2(:,1:150));
                            %                     XX = smoothdata(X,2,'gaussian',50);
                            %                                         X2 = imbilatfilt(XX,(spont_std^2)*2,2);
                            tuning{db_ind1} = X;
                            
                            X2 = imgaussfilt(X,[3,20],'Padding',0); %,'Padding','circular');
                            X2 = X2-mean2(X2(:,1:150)); % removing spont rate
                            [BF_ind, pos, ~]  = find(X2 == max(max(X2(:,200:end))));
                            
                            % checking whether peak is significant
                            %                     if mean2(X(:,200:end)) < 1 % average firing rate after onset is below 1
                            %                         BF_ind = [];
                            %                     if mean2(abs(X(:,200:end)-mean2(X(:,1:150))))<1
                            %                     BF_ind_keep = [];
                            
                            
                            BF_ind_keep = [];
                            for bff = 1:length(BF_ind)
                                if X2(BF_ind(bff),pos(bff)) > std2(X2(:,1:150))*2 % change here for BF detection threshold
                                    BF_ind_keep = [BF_ind_keep,bff];
                                end
                            end
                            BF_ind = BF_ind(BF_ind_keep);
                            pos = pos(BF_ind_keep);
                            if pos < 203 %if peak is found during the first 2ms of onset
                                BF_ind = [];
                            end
                            
                            
                            %                     if X(BF_ind,pos)<=0 && X(BF_ind,pos-1)<=0 && X(BF_ind,pos+1)<=0 % if there are no spikes within 1ms of peak latency
                            %
                            %                         BF_ind = [];
                            %                     end
                            
                            
                            %|| Calculating min latency
                            threshX3 = mean2(X(:,1:150))+std2(X(:,1:150))*2;
                            min_lat_pool = [];
                            
                            if ~isempty(BF_ind)
                                if BF_ind < 3
                                    X3 = X(BF_ind:BF_ind+2,200:end);
                                elseif BF_ind > size(X,2)-3
                                    X3 = X(BF_ind-2:BF_ind,200:end);
                                else
                                    X3 = X(BF_ind-2:BF_ind+2,200:end);
                                end
                                
                                % thresholding method 1
                                
                                X3 = mean(X3,1); % thresholding method 2 Recanzone et al 2000. comment to not use
                                
                                
                                for nx = 1:size(X3,1)
                                    ms = 1;
                                    tt = 1; %counter to break out of while loop
                                    % 3 consecutive 2ms bins must have
                                    % spikes that are over the threshold
                                    while (mean(X3(nx,ms:ms+1))<= threshX3 || ...
                                            mean(X3(nx,ms+2:ms+3))<= threshX3 || ...
                                            mean(X3(nx,ms+4:ms+5))<= threshX3) && tt~=0
                                        ms = ms+1;
                                        if ms > size(X3,2)-6
                                            
                                            tt = 0;
                                            
                                        end
                                    end
                                    if tt ~=0
                                        min_lat_pool = [min_lat_pool, ms+200];
                                    end
                                end
                                
                            end
                            
                            % Calculating min latency, end ||
                            
                            % Calculating peak firing rate with a 20ms
                            % window centered around BF peak latency.
                            X4 = [];
                            for ppp = 1:length(pos)
                                X4 = [X4 mean(X(BF_ind,pos-10:pos+10))];           
                            end
                            % Calculating peak firing rate, end
                            
                            
                            
                            if ~isempty(min_lat_pool)
                                min_lat{db_ind1} = min(min_lat_pool);
                            end
                            if figure_on == 1
                                subplot(1,4,db_ind1)
                                
                                % creating figures
                                imagesc(X2);
                                hold on
                                rectangle('Position',[200 0 100 size(X2,1)],'LineWidth',2)
                                if ~isempty(BF_ind)
                                    scatter(pos,BF_ind,100,'rp','filled')
                                    if ~isempty(min_lat_pool)
                                        scatter(min_lat{db_ind1},BF_ind,100,'cp','filled')
                                    end
                                end
                                
                                caxis([min(min(X2)),max(max(X2))])
                                title([out.tags{db_ind1+1} 'db'])
                            end
                            % creating figures end
                            
                            
                            
                            %                     BF_ind = ana_PT{db_ind}.ana.onset.BF_ind;
                            %                     FR  = ana_PT{db_ind}.ana.onset.data;
                            if ~isempty(BF_ind) && max(X4) ~=0
                                peak_FR{db_ind1} = max(X4);
                                peak_lat{db_ind1} = pos;
                                spont_FR{db_ind1} = mean2(X(:,1:150));
                                db_count = 1;
                                
                                for b = 1:length(BF_ind)
                                    out_vec{db_ind1} = [out_vec{db_ind1}; stim_freq(BF_ind(b))]; %, mean(FR(BF_ind(b),:))]];
                                    %                         out_vec = [out_vec; [nid, stim_freq(BF_ind(b)),mean(FR(BF_ind(b),:)),H_nb,tracks(t)]];
                                end
                                
                            end
                            out_stim{db_ind1} = stim_freq;
                            
                        end
                    catch
                        %                 test = 1;
                    end
                end
                if figure_on == 1
                    drawnow()
                    if db_count2 ==0
                        close(fi)
                    end
                    %             pause
                end
                
                
                if db_count ~= 0
                    db_lens = db_lens+1;
                    out.data{1}(db_lens,1) = nid;
                    out.data{1}(db_lens,2) = H_nb;
                    out.data{1}(db_lens,3) = tracks(t);
                    
                    for db_ind1 = 1:4
                        out.data{db_ind1+1}{db_lens,1} = out_vec{db_ind1};
                        out.data{db_ind1+1}{db_lens,2} = peak_lat{db_ind1};
                        out.data{db_ind1+1}{db_lens,3} = peak_FR{db_ind1};
                        out.data{db_ind1+1}{db_lens,4} = spont_FR{db_ind1};
                        out.data{db_ind1+1}{db_lens,5} = min_lat{db_ind1};
                        out.data{db_ind1+1}{db_lens,6} = tuning{db_ind1};
                        out.data{db_ind1+1}{db_lens,7} = out_stim{db_ind1};
                        
                        
                    end

                end
                
            catch
            end
            
        end
        
        %         t_ind = find(out.data{6}(:,2) == tracks(t));
        %         BF_pool = out(t_ind,2);
        %     BF_tracks(t,:) = [tracks(t),median(BF_pool)];
        % %     histogram(BF_pool,edges);
        %     hold on
        %     [N_count, edges] =histcounts(BF_pool,edges);
        %     plot(new_edges,N_count);
    end
    
end

fprintf(['done.     time : %6.2f sec \n'],toc')
% xticks([1e3 : 2.5*1e3:3*1e4])
% % set(gca,'xticklabel',num2str(get(gca,'xtick')/1e3','%.1f'))
% xticklabels(num2str(get(gca,'xtick')'/1e3))
% xlabel('kHz')
% set(gca,'xscale','log')
% % legend(num2str(tracks))
% figure
% histogram(out(:,2),edges);
% xticklabels(num2str(get(gca,'xtick')'/1e3))
% xlabel('kHz')
% set(gca,'xscale','log')



%%
% cross validation across trials
%                     total_reps = size(ana_PT{db_ind}.PSTH{1},1);
%                     BF_ind_tot = zeros(total_reps,1);
%                     pos_tot = zeros(total_reps,1);
%                     for r = 1:total_reps
%                         X = [];
%                         for st = 1:41
%                             %                         X = [X; mean(ana_PT{db_ind}.PSTH{st},1)-mean(ana_PT{db_ind}.spont)];
%                             if r ==1
%                                 X = [X; mean(ana_PT{db_ind}.PSTH{st}(r+1:end,:),1)]; %-mean(ana_PT{db_ind}.spont)]
%                             else
%                                 X = [X; mean(ana_PT{db_ind}.PSTH{st}([1:r-1,r+1:end],:),1)]; %-mean(ana_PT{db_ind}.spont)]
%                             end
%                         end
%                         XX = smoothdata(X,2,'gaussian',50);
%
%                         X2 = imgaussfilt(XX,2);
%                         X2 = X2-mean2(X2(:,1:150)); % removing spont rate
%                         [BF_ind, pos, ~]  = find(X2 == max(max(X2(:,200:end))));
%                         BF_ind_tot(r) = BF_ind;
%                         pos_tot(r) = pos;
%                     end
%                     BF_ind = floor(mean(BF_ind_tot));
%                     pos = floor(mean(pos_tot));










