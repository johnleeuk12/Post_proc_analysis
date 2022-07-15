function SUrate = get_tuning(Pool, rate, raster, figure_on)
%{

05/10/22
Despite what the function suggests, this is more than tuning:
get_tuning can plot and compare unit responses across different stimuli. 

Input:
    Pool        : pool structure from pool_data.m users should run pool_data.m as
                    below: 
                [Pool{1}, rate{1}, raster{1}] = pool_data(1,[2:3],60,1180);
                [Pool{2}, rate{2}, raster{2}] = pool_data(2320,[2:3],40);
    figure_on   : set as 1 for figures. otherwise set as 0 

Output:
    mean        : average firing rate during stimulus
    error =     : standard deviation during stim
    spont = 	: spontaneous rate
    raw =       : raw data during stimulus
    PSTH =      : raw PSTH (probably the most useful for detailed analysis
    nid =       : neuron id (to track neurons across stimuli sets
    pid =       : id corresponding to the variable Pool
    xb =        : xbz file 



%}




N_pool = size(Pool,2);
for ss = 1:N_pool
    N_list{ss} = unique([Pool{ss}.neuron_nb]);
end

% find intersecting list 
N_list_int = N_list{1};
if size(Pool,2) >1
    for ss = 1:N_pool-1
        N_list_int = intersect(N_list_int,N_list{ss+1});
    end
end
N_list_int = N_list_int(2:end);

SUrate = {};
for n = 1:length(N_list_int) % n is each neuron
    rec_list = {};
    check_counter = 1;
    SUrate{n} = {};
    
    for ss = 1:N_pool %ss is each stim for same neuron
        rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list_int(n));
        if isempty(rec_list{ss})
            check_counter = 0;
%         elseif length(rec_list{ss})<2 %This is for Puretones
%             check_counter = 0;
        end
    end
    

        
    
    if check_counter ==1
        for ss = 1:N_pool
            if ss ~=3
            p = rec_list{ss}(end); % try end as well change for figures
            else
                p = rec_list{ss}(end);
            end
            nreps = size(rate{ss}.stim{p},2);
            SUrate{n}{ss}.mean = mean(rate{ss}.stim{p},2);
            SUrate{n}{ss}.error = std(rate{ss}.stim{p},1,2)/sqrt(nreps);
            SUrate{n}{ss}.spont = mean(rate{ss}.pre{p},2);
            SUrate{n}{ss}.raw = rate{ss}.stim{p};
            SUrate{n}{ss}.PSTH = {};
            SUrate{n}{ss}.nid = Pool{ss}(p).neuron_nb;
            SUrate{n}{ss}.pid = p;
            SUrate{n}{ss}.xb = Pool{ss}(p).xb;
            SUrate{n}{ss}.xb.data = [];
            
            for te = 1:size(rate{ss}.PSTH,2)
                SUrate{n}{ss}.PSTH{te}= rate{ss}.PSTH{p,te};
            end
            StimDur = Pool{ss}(p).xb.stimulus_ch1(:,5)*1e-3;
            
            % Change stim label
            if Pool{ss}(p).xb.analysis_code >100 && Pool{ss}(p).xb.analysis_code< 1000
                plot_type{ss} = 'User';
                stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,1);
                for st = 1: length(stim_label{ss}{p})
                    idx1 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},': ')+2;
                    idx2 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
                    idx3 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'rev');
                    idx4 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'.txt');
                    idx5 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'_t');
                    
                    %                     stim_name{ss}{st} = Pool{ss}(p).xb.user_stimulus_desc_ch1{st}(idx1:idx4);
                    stim_name{ss}{st} = [Pool{ss}(p).xb.user_stimulus_desc_ch1{st}(idx1:idx5) Pool{ss}(p).xb.user_stimulus_desc_ch1{st}(idx3:idx4)];
                    
                    
                end
            elseif  strcmp(Pool{ss}(p).xb.analysis_type,'User: Trill_list.txt')
                plot_type{ss} = 'User';
                stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,1);
                for st = 1: length(stim_label_1)
                    idx1 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},': ')+2;
                    idx2 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
                    idx3 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'rev');
                    idx5 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'_t');
                    idx4 = strfind(Pool{ss}(p).xb.user_stimulus_desc_ch1{st},'.txt');
                    stim_name{ss}{st} = [Pool{ss}(p).xb.user_stimulus_desc_ch1{st}(idx1:idx5) Pool{ss}(p).xb.user_stimulus_desc_ch1{st}(idx3:idx4)];
                    %             if isempty(idx3) ~=1
                    %                 stim_name{st} = [stim_name{st} '_rev'];
                    %             end
                    
                end
                
                %     elseif Pool_1(p).xb.analysis_code >2000
                
            elseif      Pool{ss}(p).xb.analysis_code == 1
                plot_type{ss} = 'Tones';
                stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
                max_id = find(SUrate{n}{ss}.mean == max(SUrate{n}{ss}.mean));
                
%                 %Calculating Best frequency for tones
%                 if SUrate{n}{ss}.mean(max_id(1)) > SUrate{n}{ss}.spont +2*SUrate_1{n}.error(max_id(1))
%                     BF_pool = [BF_pool stim_label_1(max_id(1))];
%                     Hole_Track = [Hole_Track; Pool_1(p).hole_nb Pool_1(p).track_nb];
%                 end

            elseif Pool{ss}(p).xb.analysis_code == 62
                plot_type{ss} = 'Noise';
                stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
                max_id = find(SUrate{n}{ss}.mean == max(SUrate{n}{ss}.mean));
            else
                plot_type{ss} = 'VT';
                %                 stim_label{ss} = Pool{ss}(p).xb.stimulus_ch1(:,9); %change, modify, make better
                [stim_label{ss}{p}, stim_ticks{ss}{p}] = stim_index(Pool{ss},p);
                
            end
            %
            nreps = Pool{ss}(p).xb.stimulus_ch1(1,4);
            nStim = max(Pool{ss}(p).xb.stimulus_ch1(:,1));

            TotalReps = nStim*nreps;
            
            % for real vocalizations
            %     nat_rev(n,1) = mean(SUrate_1{n}.mean(1:2:20)); %natural
            %     nat_rev(n,2) = mean(SUrate_1{n}.mean(2:2:20)); %reversed
            
   
        end
        
        
        %Figures
        if figure_on == 1
            pause(0.1)
            %                 if SUrate{n}{ss}.spont >0 &&
            figure(n)
            set(gcf, 'Position', [400 100 1700 800]);
            for ss = 1:N_pool
                p = rec_list{ss}(1); % try end as well
                
                %Tuning curve
                subplot(2,N_pool,ss)
                errorbar(stim_label{ss}{p},SUrate{n}{ss}.mean,SUrate{n}{ss}.error,'LineWidth',2);
                hold on
                %             plot(stim_label{ss},ones(1,length(stim_label{ss}))*SUrate{n}{ss}.spont,'--k');
                plot(stim_label{ss}{p},ones(1,length(stim_label{ss}{p}))*mean(SUrate{n}{ss}.spont),'--k')
                if strcmp(plot_type{ss},'User')
                    xticks(stim_label{ss}{p});
                    xticklabels(stim_name{ss})
                    xtickangle(45)
                elseif strcmp(plot_type{ss},'Tones')
                    set(gca,'xscale','log')
                    xticks(round((stim_label{ss}{p}(1:5:end).')*1e-1)*1e1);
                    xticklabels(round(stim_label{ss}{p}(1:5:end)*1e-1)*1e-2)
                    
                    xtickangle(45)
                    xlabel('kHz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                    
                elseif strcmp(plot_type{ss},'Noise')
                    set(gca,'xscale','log')
                    xticks(round((stim_label{ss}{p}(1:5:end).')*1e1)*1e-1);
                    %                 xticklabels(round(stim_label{ss}(1:5:end)*1e1)*1e-1)
                    xtickangle(45)
                    xlabel('kHz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                end
                
                
                
                %rasterplot
                subplot(2,N_pool,ss+N_pool)
                
                PreStim = Pool{ss}(p).xb.pre_stimulus_record_time*1e-3; %s
                PostStim = Pool{ss}(p).xb.post_stimulus_record_time*1e-3; %s
                StimDur = Pool{ss}(p).xb.stimulus_ch1(:,5)*1e-3;
                
                nreps = Pool{ss}(p).xb.stimulus_ch1(1,4);
                nStim = max(Pool{ss}(p).xb.stimulus_ch1(:,1));
                
                TotalReps = nStim*nreps;
                
                
                
                for st = 1:length(StimDur)
                    rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                end
                hold on
                plot(raster{ss}.spikes{p},nreps*(raster{ss}.stim{p}-1)+raster{ss}.rep{p},'k.','MarkerSize',15);
                %     pause
                xlabel('time (s)')
                
                yticks([1:nreps*2:TotalReps]+10)
                %             yticks([1:nreps*2:TotalReps]+10)
                %             stim_ticks{ss} = {};
                
                %             for stim = 1:length(stim_label{ss})
                %                 %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                %                 stim_ticks{ss}{stim}=num2str(round(stim_label{ss}(stim)*10)/10);
                %             end
                if strcmp(plot_type{ss},'User')
                    
                    yticks([1:nreps:TotalReps]+5)
                    
                    yticklabels(stim_name{ss})
                elseif strcmp(plot_type{ss},'VT')
                    ylabel('kHz')
                    yticklabels(stim_ticks{ss}{p}(1:2:end))
                else
                    stim_ticks{ss} = {};
                    
                    for stim = 1:length(stim_label{ss}{p})
                        %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                        stim_ticks{ss}{p}{stim}=num2str(round(stim_label{ss}{p}(stim)*10)/10);
                    end
                    yticklabels(stim_ticks{ss}{p}(1:2:end))
                end
                axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                hold off
                title('rasterplot')
                
                
                
                
                drawnow
            end
            sgtitle(['Ch' num2str(Pool{ss}(p).best_ch)...
                ' H' num2str(Pool{ss}(p).hole_nb) 'T' ...
                num2str(Pool{ss}(p).track_nb) ' Unit ' num2str(Pool{ss}(p).neuron_nb)])
        end
    end
end



 