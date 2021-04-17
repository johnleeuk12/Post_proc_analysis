function out = ana_BF(H_nb, figure_on)

%% CF analysis 03/07/2021
% load('M60F_unit_list_new.mat');
animal_name = 'M60f';
% H_nb = 4;
% db_ind = 2;4 is 20db % 3 is 40db, 2 is 60db 1 is 80 db

edges = logspace(3,4.5,25);
% new_edges = zeros(1,39);
% for e = 1:39
%     new_edges(e) = (edges(e) + edges(e+1))/2;
% end



savedir = fullfile('E:\DATA\ana_tones',filesep,animal_name);
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
hole_ind = find([neurons_loc_tag.hole_nb] == H_nb);
tracks = unique([neurons_loc_tag(hole_ind).track_nb]);
N_list = [neurons_loc_tag.neuron_nb];
BF_tracks = [];
% figure
db_lens = 0;
for t = 1:length(tracks)
    sub_n_list = N_list(find([neurons_loc_tag.track_nb] ==tracks(t) &[neurons_loc_tag.hole_nb] == H_nb ));
    for n =1:length(sub_n_list)
        nid = sub_n_list(n);
        PTname= 'PTn0000';
        PTname = [PTname(1:end-length(num2str(nid))) num2str(nid) '.mat'];
        load(PTname);
        ana_PT = filename; %currently loaded file is filename. change afterwards
        db_count = 0; % counter to see if there's at least one non-empty BF_ind
        out_vec = {};
        % figures
        if figure_on ==  1
            fi = figure;
            set(fi, 'Position', [100 200 2400 660]);
            sgtitle(['unit' num2str(nid) ' H4' 'T' num2str(tracks(t))])
        end
        % figures end
        db_count2 = 0;
        %         db_count3 = 0;
        for db_ind = 1:4
            out_vec{db_ind} = [];
            try
                if ~isempty(ana_PT{db_ind})
                    db_count2 = 1;
                    

                    
                    X = [];
                    for st = 1:41
                        %                         X = [X; mean(ana_PT{db_ind}.PSTH{st},1)-mean(ana_PT{db_ind}.spont)];
                        X = [X; mean(ana_PT{db_ind}.PSTH{st},1)]; %-mean(ana_PT{db_ind}.spont)]
                        
                    end
                    
                    %                     spont_std = std2(X2(:,1:150));
                    XX = smoothdata(X,2,'gaussian',50);
                    %                     X2 = imbilatfilt(XX,(spont_std^2)*2,2);
                    X2 = imgaussfilt(XX,2);
                    X2 = X2-mean2(X2(:,1:150)); % removing spont rate
                    [BF_ind, pos, ~]  = find(X2 == max(max(X2(:,200:end))));

                    % checking whether peak is significant
                    if mean2(X(:,200:end)) < 1 % average firing rate is below 0.5
                        BF_ind = [];
%                     if mean2(abs(X(:,200:end)-mean2(X(:,1:150))))<1
                    elseif X2(BF_ind,pos) < std2(X2(:,1:150))*2
                        BF_ind = [];
                    end
                    
                    
                    if figure_on == 1
                        subplot(1,4,db_ind)
                        
                        % creating figures
                        imagesc(X2);
                        hold on
                        rectangle('Position',[200 0 100 size(X2,1)],'LineWidth',2)
                        if ~isempty(BF_ind)
                            scatter(pos,BF_ind,100,'rp','filled')
                        end
                        
                        caxis([min(min(X2)),max(max(X2))])
                        title([out.tags{db_ind+1} 'db'])
                    end
                    % creating figures end
                    
                    
                    
                    %                     BF_ind = ana_PT{db_ind}.ana.onset.BF_ind;
                    %                     FR  = ana_PT{db_ind}.ana.onset.data;
                    if ~isempty(BF_ind)
                        db_count = 1;
                        
                        stim_freq = ana_PT{db_ind}.xb.stimulus_ch1(:,8);
                        for b = 1:length(BF_ind)
                            out_vec{db_ind} = [out_vec{db_ind}; stim_freq(BF_ind(b))]; %, mean(FR(BF_ind(b),:))]];
                            %                         out_vec = [out_vec; [nid, stim_freq(BF_ind(b)),mean(FR(BF_ind(b),:)),H_nb,tracks(t)]];
                        end
                    end
                end
            catch
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
            for db_ind = 1:4
                out.data{db_ind+1}{db_lens,1} = out_vec{db_ind};
            end
            out.data{6}(db_lens,:) = [H_nb, tracks(t)];
        end
        
        
%       test=  1;
      
    end
    
    %         t_ind = find(out.data{6}(:,2) == tracks(t));
    %         BF_pool = out(t_ind,2);
    %     BF_tracks(t,:) = [tracks(t),median(BF_pool)];
    % %     histogram(BF_pool,edges);
    %     hold on
    %     [N_count, edges] =histcounts(BF_pool,edges);
    %     plot(new_edges,N_count);
end
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










