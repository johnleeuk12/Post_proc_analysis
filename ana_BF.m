function [out, BF_tracks] = ana_BF(H_nb,db_ind)

%% CF analysis 03/07/2021
% load('M60F_unit_list_new.mat');
animal_name = 'M60f';
H_nb = 4;
db_ind = 2; % 3 is 40db, 2 is 60db 1 is 80 db

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


out = [];
hole_ind = find([neurons_loc_tag.hole_nb] == H_nb);
tracks = unique([neurons_loc_tag(hole_ind).track_nb]);
N_list = [neurons_loc_tag.neuron_nb];
BF_tracks = [];
% figure
for t = 1:length(tracks)
    sub_n_list = N_list(find([neurons_loc_tag.track_nb] ==tracks(t) &[neurons_loc_tag.hole_nb] == H_nb ));
    for n =1:length(sub_n_list)
        nid = sub_n_list(n);
        PTname= 'PTn0000';
        PTname = [PTname(1:end-length(num2str(nid))) num2str(nid) '.mat'];
        load(PTname);
        ana_PT = filename; %currently loaded file is filename. change afterwards
        if ~isempty(ana_PT{db_ind})
            BF_ind = ana_PT{db_ind}.ana.onset.BF_ind;
            FR  = ana_PT{db_ind}.ana.onset.data;
            if ~isempty(BF_ind)
                stim_freq = ana_PT{db_ind}.xb.stimulus_ch1(:,8);
                for b = length(BF_ind)
                    out = [out; [nid, stim_freq(BF_ind(b)),mean(FR(BF_ind(b),:)),H_nb,tracks(t)]];
                end
            end
        end
    end
    t_ind = find(out(:,5) == tracks(t));
    BF_pool = out(t_ind,2);
    BF_tracks(t,:) = [tracks(t),median(BF_pool)];
%     histogram(BF_pool,edges);
    hold on
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
