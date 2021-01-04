
Pool_1 = Pool{1};
raster_1 = raster{1};
rate_1 = rate{1};

% for n = 79:88
% Pool_1(n).track_nb = 1;
% end
figure_on = 0;
N_list_1 = unique([Pool_1.neuron_nb]);
% N_list_2 = unique([Pool_2.neuron_nb]);
N_list_1 = N_list_1(2:end);
SUrate_1 = {};
BF_pool = [];
Latency_pool = [];
Hole_Track = [];
nat_rev = [];
for n = 1:length(N_list_1)
    rec_list_1 = find([Pool_1.neuron_nb] ==N_list_1(n));
    p = rec_list_1(1); %try end as well
    nreps = size(rate_1.stim{1},2);
    SUrate_1{n}.mean = mean(rate_1.stim{p},2); % -mean(rate_1.pre{p},2); %(Spikes/second)
    SUrate_1{n}.error = std(rate_1.stim{p},1,2)/sqrt(nreps);
    SUrate_1{n}.spont = mean(mean(rate_1.pre{p}));
    SUrate_1{n}.raw = rate_1.stim{p};
    StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;
    
    % Change stim label
    if Pool_1(p).xb.analysis_code >99 && Pool_1(p).xb.analysis_code< 1000
        plot_type = 'User';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,1);
        for st = 1: length(stim_label_1)
            idx1 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},': ')+2;
            idx2 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
            idx3 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'rev');
            idx4 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'.txt');
            stim_name{st} = Pool_1(p).xb.user_stimulus_desc_ch1{st}(idx1:idx4);
%             if isempty(idx3) ~=1
%                 stim_name{st} = [stim_name{st} '_rev'];
%             end
            
        end
    elseif  strcmp(Pool_1(p).xb.analysis_type,'User: Trill_list.txt') 
        plot_type = 'User';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,1);
        for st = 1: length(stim_label_1)
            idx1 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},': ')+2;
            idx2 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
            idx3 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'rev');
            idx4 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'.txt');
            stim_name{st} = Pool_1(p).xb.user_stimulus_desc_ch1{st}(idx1:idx4);
%             if isempty(idx3) ~=1
%                 stim_name{st} = [stim_name{st} '_rev'];
%             end
            
        end
         
        %     elseif Pool_1(p).xb.analysis_code >2000
        
    elseif      Pool_1(p).xb.analysis_code == 1
        plot_type = 'Tones';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,8);
        max_id = find(SUrate_1{n}.mean == max(SUrate_1{n}.mean));
        
        %Calculating Best frequency for tones
        if SUrate_1{n}.mean(max_id(1)) > SUrate_1{n}.spont +2*SUrate_1{n}.error(max_id(1))
            BF_pool = [BF_pool stim_label_1(max_id(1))];
            
            Hole_Track = [Hole_Track; Pool_1(p).hole_nb Pool_1(p).track_nb];
            %Calculating Stim Latency.
            PSTH = mean(rate_1.PSTH{p,max_id(1)},1);
            latency = find(PSTH(Pool_1(p).xb.pre_stimulus_record_time:end)> 2*std2(rate_1.pre{p}));
            latency = latency(1); %first crossing
            Latency_pool = [Latency_pool, latency];
        end
        

    else
        plot_type = 'VT';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,9); %change, modify, make better
        
    end
    %
    nreps = Pool_1(p).xb.stimulus_ch1(1,4);
    nStim = max(Pool_1(p).xb.stimulus_ch1(:,1));
    %             nreps = 10;
    
    %
    % nStim = nStim + max(x2.stimulus_ch1(:,1));
    %
    TotalReps = nStim*nreps;
    
    % for real vocalizations
%     nat_rev(n,1) = mean(SUrate_1{n}.mean(1:2:20)); %natural
%     nat_rev(n,2) = mean(SUrate_1{n}.mean(2:2:20)); %reversed
%     
    
    if figure_on == 1
        if SUrate_1{n}.spont >0% && SUrate_2{n}.spont >1
            
            
            
            figure
            set(gcf, 'Position', [400 100 1700 800]);
            
            
            subplot(1,2,1)
            errorbar(stim_label_1,SUrate_1{n}.mean,SUrate_1{n}.error,'LineWidth',2);
            hold on
            plot(stim_label_1,ones(1,length(stim_label_1))*SUrate_1{n}.spont,'--k');
            
            if strcmp(plot_type,'User')
                xticks(stim_label_1);
                xticklabels(stim_name)
                xtickangle(45)
            elseif strcmp(plot_type,'Tones')
                set(gca,'xscale','log')
                xticks(round((stim_label_1(1:4:end).')*1e-1)*1e1);
                xticklabels(round(stim_label_1(1:4:end)*1e-1)*1e-2)
                
                xtickangle(45)
                xlabel('kHz')
                ylabel('firing rate (spikes/s)')
                title('tuning curve')
                
            end
            
            
            
            subplot(1,2,2)
            PreStim = Pool_1(p).xb.pre_stimulus_record_time*1e-3; %s
            PostStim = Pool_1(p).xb.post_stimulus_record_time*1e-3; %s
            StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;
            
            nreps = Pool_1(p).xb.stimulus_ch1(1,4);
            nStim = max(Pool_1(p).xb.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            
            
            hold off
            for st = 1:length(StimDur)
                rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold on
            plot(raster_1.spikes{p},nreps*(raster_1.stim{p}-1)+raster_1.rep{p},'k.','MarkerSize',15);
            hold on
            
            if ~isempty(Latency_pool)
                if ~isempty(latency)
                    scatter(latency*1e-3,max_id(1)*nreps,300,'rx');
                end
            end
            %     pause
            xlabel('time (s)')
            ylabel('Hz')
            yticks([1:nreps*2:TotalReps]+10)
            %             yticks([1:nreps*2:TotalReps]+10)
            stim_ticks = {};
            
            for stim = 1:length(stim_label_1)
                %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                stim_ticks{stim}=num2str(round(stim_label_1(stim)*10)/10);
            end
            if strcmp(plot_type,'User')
                
                yticks([1:nreps:TotalReps]+5)
                
                yticklabels(stim_name)
            else
                yticklabels(stim_ticks(1:2:end))
            end
            axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
            hold off
            title('rasterplot')
            
            sgtitle(['Ch' num2str(Pool_1(p).best_ch) ' H' num2str(Pool_1(p).hole_nb) 'T' num2str(Pool_1(p).track_nb) ' u' num2str(Pool_1(p).neuron_nb)])
            
            
            drawnow
        end
    end
    
end

% %% select twitter responsive sessions 
% 
% % not all neurons respond to all twitters 
% stim_length = length(SUrate_1{1}.mean);
% sig_list = zeros(1,length(N_list_1));
% for n = 1:length(N_list_1)
%     if sum(SUrate_1{n}.mean>SUrate_1{n}.spont+ 2*SUrate_1{n}.error)
%         sig_list(n) = 1;
%     end
% end
% 
% sig_list_id = find(sig_list == 1);
% 
% % sig_list_id = 1:48;
% 
% 
% %IPI twitter study
% pop_raw = [];
% figure
% for i = 1:length(sig_list_id)
%     
% %     errorbar(StimDur,SUrate_1{sig_list_id(i)}.mean, SUrate_1{sig_list_id(i)}.error)
%     hold on
%     pop_raw = [pop_raw SUrate_1{sig_list_id(i)}.mean/max(SUrate_1{sig_list_id(i)}.mean)];
% end
% 
% figure
% pop_ave = mean(pop_raw,2);
% pop_error = std(pop_raw,[],2)/sqrt(length(sig_list_id));
% 
% errorbar(StimDur,mean(pop_raw,2),pop_error)
% 
% 
% figure
% 
% scatter(nat_rev(:,1),nat_rev(:,2))
% hold on
% plot([0:80],[0:80])
% 
% test = nat_rev(:,1)./nat_rev(:,2);
% edges = [0.5:0.1:1.5];
% histogram(test,edges);
% 
% 
% 
% 
% 
% natVSrev = [];
% for p = 1:length(sig_list_id)
%     n = sig_list_id(p);
%     for m = 1:round(stim_length/2)
%     natVSrev(m,p) = SUrate_1{n}.mean(2*m-1)/SUrate_1{n}.mean(2*m);
%     if SUrate_1{n}.mean(2*m) ==  0 && SUrate_1{n}.mean(2*m-1) ==  0
%         natVSrev(m,p) = 1;
%     elseif SUrate_1{n}.mean(2*m) ==  0 && SUrate_1{n}.mean(2*m-1) ~=  0
%         natVSrev(m,p) = 0.01;
%     elseif SUrate_1{n}.mean(2*m-1) ==  0 && SUrate_1{n}.mean(2*m) ~=  0
%          natVSrev(m,p) = 10;
%     end 
%     end
% end
% 
% histogram(natVSrev)
% 
% test = mean(natVSrev,1);
% 
% 
% histogram(log10(test))
% edges = [-1.3:0.1:1];
% histogram(log10(test),edges)
% hold on
% % histogram(log10(test_2),edges)
% 
% % test = nat_rev(:,1)./nat_rev(:,2); 
% 
% 
% % 
% % test = log10(test);
% % test = rmmissing(test); % above matlab 2018b
% % 
% % [P,h,stats] = signrank(test)
% % 
% % med_test = median(test);
% % 
% 
% % 
% % hold on
% % % plot([1:20],[1:20])
% edges = logspace(2,5,31);
% % % % histogram(stim_label(BF_pool),edges);
% histogram(BF_pool,edges) %,stim_label_1)
% set(gca,'xscale','log')
% xticklabels({0.1, 1, 10})
% xlabel('khz')

