function det_BF(Pool, rate, raster)


%%




animal_name = 'M60F';
savedir = fullfile('E:\DATA\ana_tones',filesep,animal_name);

try 
    load([savedir '\neurons_loc_tag.mat']);
catch
    neurons_loc_tag = {};
    save([savedir '\neurons_loc_tag.mat'], 'neurons_loc_tag');
end





N_list = unique([Pool.neuron_nb]);
N_list = N_list(2:end);
SUrate = {};

for n = 1:length(N_list)
        rec_list = find([Pool.neuron_nb] == N_list(n));
        p = rec_list(1);
        neurons_loc_tag(Pool(p).neuron_nb).neuron_nb = Pool(p).neuron_nb;
        neurons_loc_tag(Pool(p).neuron_nb).hole_nb = Pool(p).hole_nb;
        neurons_loc_tag(Pool(p).neuron_nb).track_nb = Pool(p).track_nb;
end

save([savedir '\neurons_loc_tag.mat'], 'neurons_loc_tag');

T_list = unique([neurons_loc_tag.neuron_nb]);


tic
for n = 1:length(N_list)
    rec_list = find([Pool.neuron_nb] == N_list(n));
    check_counter = 1;
    SUrate{n} = {};
    
    for r = 1: length(rec_list)
        p= rec_list(r);
        a_code = Pool(p).xb.analysis_code;
        db = Pool(p).xb.stimulus_ch1(1,3);
        lens = Pool(p).xb.stimulus_ch1(1,5);
        if lens ==100
            if a_code == 1
                switch db
                    case 80
                        s2 = 1;
                    case 60
                        s2 = 2;
                    case 40
                        s2 = 3;
                    case 20
                        s2 = 4;
                end
            elseif a_code == 62
                switch db
                    case 80
                        s2 = 5;
                    case 60
                        s2 = 6;
                    case 40
                        s2 = 7;
                    case 20
                        s2 = 8;
                end
            end
            s1 = 2;
            try if isempty(SUrate{n}{1,s2})
                    s1 = 1;
                end
            catch
                s1 = 1;
            end
            
            nreps = size(rate.stim{p},2);
            SUrate{n}{s1,s2}.mean = mean(rate.stim{p},2);
            %                 SUrate{n}{s1,s2}.post = mean(rate.post{p},2);
            %                 SUrate{n}{s1,s2}.error = std(rate.stim{p},1,2)/sqrt(nreps);
            SUrate{n}{s1,s2}.spont = mean(rate.pre{p},2);
            %                 SUrate{n}{s1,s2}.raw = rate.stim{p};
            SUrate{n}{s1,s2}.PSTH = {};
            SUrate{n}{s1,s2}.nid = Pool(p).neuron_nb;
            SUrate{n}{s1,s2}.pid = p;
            SUrate{n}{s1,s2}.xb = Pool(p).xb;
            SUrate{n}{s1,s2}.xb.data = [];
            
            for te = 1:length(SUrate{n}{s1,s2}.mean)
                SUrate{n}{s1,s2}.PSTH{te}= rate.PSTH{p,te};
            end
            
            % Determine peaks and latency, only for s1 = 1;
            SUrate{n}{s1,s2}.ana = {};
            SUrate{n}{s1,s2}.ana.onset.data = zeros(length(SUrate{n}{s1,s2}.PSTH),nreps);
            SUrate{n}{s1,s2}.ana.offset1.data = zeros(length(SUrate{n}{s1,s2}.PSTH),nreps);
            SUrate{n}{s1,s2}.ana.offset2.data = zeros(length(SUrate{n}{s1,s2}.PSTH),nreps);
            
            for st = 1:length(SUrate{n}{s1,s2}.PSTH)
                SUrate{n}{s1,s2}.ana.onset.data(st,:) = mean(SUrate{n}{s1,s2}.PSTH{st}(:,200:300),2).';
                SUrate{n}{s1,s2}.ana.offset1.data(st,:) = mean(SUrate{n}{s1,s2}.PSTH{st}(:,300:400),2).';
                SUrate{n}{s1,s2}.ana.offset2.data(st,:) = mean(SUrate{n}{s1,s2}.PSTH{st}(:,400:500),2).';
%                 SUrate{n}{s1,s2}.ana.spont.data = 
            end
            
            % onset
            temp_mean = movmean(mean(SUrate{n}{s1,s2}.ana.onset.data,2),3);            
            [pks,loc]= findpeaks(temp_mean);
            bf_ind = loc(find(pks ==max(pks)));
            SUrate{n}{s1,s2}.ana.onset.BF_ind = [];
            for s_bf = 1:length(bf_ind)
                if temp_mean(bf_ind(s_bf))> mean(SUrate{n}{s1,s2}.spont,1) + 2*std(SUrate{n}{s1,s2}.spont,[],1) && temp_mean(bf_ind(s_bf))>2
                    SUrate{n}{s1,s2}.ana.onset.BF_ind = [SUrate{n}{s1,s2}.ana.onset.BF_ind bf_ind(s_bf)];
                    
                end
            end
            % offset1
            temp_mean = movmean(mean(SUrate{n}{s1,s2}.ana.offset1.data,2),3);
            [pks,loc]= findpeaks(temp_mean);
            bf_ind = loc(find(pks ==max(pks)));
            SUrate{n}{s1,s2}.ana.offset1.BF_ind = [];
            for s_bf = 1:length(bf_ind)
                if temp_mean(bf_ind(s_bf))> mean(SUrate{n}{s1,s2}.spont,1) + 2*std(SUrate{n}{s1,s2}.spont,[],1) && temp_mean(bf_ind(s_bf))>2
                    SUrate{n}{s1,s2}.ana.offset1.BF_ind = [SUrate{n}{s1,s2}.ana.offset1.BF_ind bf_ind(s_bf)];
                    
                end
            end
            % offset2
            SUrate{n}{s1,s2}.ana.offset2.BF_ind = [];
            temp_mean = movmean(mean(SUrate{n}{s1,s2}.ana.offset2.data,2),3);
            
            [pks,loc]= findpeaks(temp_mean);
            bf_ind = loc(find(pks ==max(pks)));
            for s_bf = 1:length(bf_ind)
                if temp_mean(bf_ind(s_bf))> mean(SUrate{n}{s1,s2}.spont,1) + 2*std(SUrate{n}{s1,s2}.spont,[],1) && temp_mean(bf_ind(s_bf))>2
                    SUrate{n}{s1,s2}.ana.offset2.BF_ind = [SUrate{n}{s1,s2}.ana.offset2.BF_ind bf_ind(s_bf)];
                    
                end
            end
            nid = SUrate{n}{s1,s2}.nid;
            
        end
    end
    PTname= 'PTn0000';
    PTname = [PTname(1:end-length(num2str(nid))) num2str(nid) '.mat'];
    filepath = fullfile(savedir, filesep, PTname);
    if mod(n,10) ==1
        fprintf(['%4d /' num2str(length(N_list)) ' time : %6.2f sec \n'],n,toc')
    end
    ana_PT = SUrate{n};
    save(filepath,'ana_PT');
    
end

%% Determine peaks and latency

% for n = 1%:length(N_list)
% 
%     
%     
% 
% end
% 



%% plotting figures

for n = 1:10 %80 % 1:length(N_list)
    fi = figure(n);
    %         set(fi, 'Position', [1 1 2559 1060]);
    fi.WindowState = 'maximized';
    for s1 =1:2
        for s2 =1:4
            try if ~isempty(SUrate{n}{s1,s2})
                    %                     if s2<4
                    subplot(2,8,(s1-1)*8+1+(s2-1)*2+1)
                    %                     else
                    %                         subplot(2,6,(s1-1)*16+1+(4-1)*2+1)
                    %                     end
                    hold off
                    stim_label = SUrate{n}{s1,s2}.xb.stimulus_ch1(:,8);
                    stim_ticks = {};
                    
                    for stim = 1:length(stim_label)
                        stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                    end
                    
                    p  =SUrate{n}{s1,s2}.pid;
                    %tuning curve
                    %                     post_stim_resp = [];
                    %                     post_stim_err = [];
                    %                     for st = 1:length(SUrate{n}{s1,s2}.PSTH)
                    %                         post_stim_resp(st) = mean2(SUrate{n}{s1,s2}.PSTH{st}(:,300:500));
                    %                         post_stim_err(st) = std2(SUrate{n}{s1,s2}.PSTH{st}(:,300:500))/sqrt(1e3);
                    %                     end
                    %                     max_ind_post = find(movmean(post_stim_resp)
                    temp_mean1 = movmean(mean(SUrate{n}{s1,s2}.ana.onset.data,2),3);
                    temp_mean2 = movmean(mean(SUrate{n}{s1,s2}.ana.offset1.data,2),3);
                    temp_mean3 = movmean(mean(SUrate{n}{s1,s2}.ana.offset2.data,2),3);
                    plot(temp_mean1,stim_label,'-r','LineWidth',2);
                    %                     errorbar(SUrate{n}{s1,s2}.mean,stim_label,SUrate{n}{s1,s2}.error/sqrt(5),'horizontal','-r','LineWidth',2);
                    
                    hold on
                    plot(temp_mean2,stim_label,'-b','LineWidth',2);
                    plot(temp_mean3,stim_label,'-g','LineWidth',2);
                    if ~isempty(SUrate{n}{s1,s2}.ana.onset.BF_ind)
                        scatter(temp_mean1(SUrate{n}{s1,s2}.ana.onset.BF_ind),stim_label(SUrate{n}{s1,s2}.ana.onset.BF_ind),'*r');
                    end
                    if ~isempty(SUrate{n}{s1,s2}.ana.offset1.BF_ind)
                        scatter(temp_mean2(SUrate{n}{s1,s2}.ana.offset1.BF_ind),stim_label(SUrate{n}{s1,s2}.ana.offset1.BF_ind),'*b');
                    end
                    if ~isempty(SUrate{n}{s1,s2}.ana.offset2.BF_ind)
                        scatter(temp_mean3(SUrate{n}{s1,s2}.ana.offset2.BF_ind),stim_label(SUrate{n}{s1,s2}.ana.offset2.BF_ind),'*g');
                    end
                    %                     plot(movmean(post_stim_resp,5),stim_label,'-b','LineWidth',2);
                    %                     errorbar(post_stim_resp,stim_label,post_stim_err,'horizontal','-b','LineWidth',2);
                    
                    xline(mean(SUrate{n}{s1,s2}.spont),'--k');
%                     legend('stim','post stim','spont')
                    set(gca, 'YScale', 'log')
                    axis([0 inf min(stim_label) max(stim_label)])
                    title(['Tuning,' num2str(SUrate{n}{s1,s2}.xb.stimulus_ch1(1,3)) 'db'])
                    % rasterplot
                    subplot(2,8,(s1-1)*8+1+(s2-1)*2)
                    hold off
                    PreStim = SUrate{n}{s1,s2}.xb.pre_stimulus_record_time*1e-3; %s
                    PostStim = SUrate{n}{s1,s2}.xb.post_stimulus_record_time*1e-3; %s
                    StimDur = SUrate{n}{s1,s2}.xb.stimulus_ch1(:,5)*1e-3;
                    
                    nreps = SUrate{n}{s1,s2}.xb.stimulus_ch1(1,4);
                    nStim = max(SUrate{n}{s1,s2}.xb.stimulus_ch1(:,1));
                    
                    TotalReps = nStim*nreps;
                    
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{p},nreps*(raster.stim{p}-1)+raster.rep{p},'k.','MarkerSize',15);
                    %     pause
                    xlabel('time (s)')
                    
                    yticks([1:nreps*2:TotalReps]+10)
                    db = SUrate{n}{s1,s2}.xb.stimulus_ch1(1,3);
                    if SUrate{n}{s1,s2}.xb.analysis_code == 1
                        title(['PT ' num2str(db) 'db att']);
                    else
                        title(['noise' num2str(db) 'db att']);
                    end
                    
                    
                    yticklabels(stim_ticks(1:2:end))
                    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                    hold off
                    title('rasterplot')
                    
                    
                end
            catch
            end
            
        end
    end
    drawnow
    
    pause(0.1);
    
    sgtitle([' H' num2str(Pool(p).hole_nb) 'T' ...
        num2str(Pool(p).track_nb) ' Unit ' num2str(Pool(p).neuron_nb)])
    
end




















