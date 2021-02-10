function SUplot(nid_list,SUrate,Pool, raster)
%%

% nid_list = [35,36,158];


nid_all = zeros(1,length(SUrate));
for n = 1:length(SUrate)
    nid_all(n) = SUrate{n}{1}.nid;
end

for ni = 1:length(nid_list)
    figure(ni)
    for ss = 1:length(SUrate{1})
        nid = nid_list(ni);
        n = find(nid_all == nid_list(ni));
        p = SUrate{n}{ss}.pid;
        
        %plotting figures
        
        StimDur = Pool{ss}(p).xb.stimulus_ch1(:,5)*1e-3;
        nreps = Pool{ss}(p).xb.stimulus_ch1(1,4);
        nStim = max(Pool{ss}(p).xb.stimulus_ch1(:,1));
        TotalReps = nStim*nreps;
        
        if      Pool{ss}(p).xb.analysis_code == 1
            plot_type{ss} = 'Tones';
            stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
            %             max_id = find(SUrate{n}{ss}.mean == max(SUrate{n}{ss}.mean));
            
            
        elseif Pool{ss}(p).xb.analysis_code == 62
            plot_type{ss} = 'Noise';
            stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
            %             max_id = find(SUrate{n}{ss}.mean == max(SUrate{n}{ss}.mean));
        else
            plot_type{ss} = 'VT';
            %                 stim_label{ss} = Pool{ss}(p).xb.stimulus_ch1(:,9); %change, modify, make better
            [stim_label{ss}{p}, stim_ticks{ss}{p}] = stim_index(Pool{ss},p);
            
        end
        
        set(gcf, 'Position', [400 100 1700 800]);

%         p = rec_list{ss}(1); % try end as well
        
        %Tuning curve
        subplot(2,length(SUrate{1}),ss)
        errorbar(stim_label{ss}{p},SUrate{n}{ss}.mean,SUrate{n}{ss}.error,'LineWidth',2);
        hold on
        %             plot(stim_label{ss},ones(1,length(stim_label{ss}))*SUrate{n}{ss}.spont,'--k');
        plot(stim_label{ss}{p},ones(1,length(stim_label{ss}{p}))*mean(SUrate{n}{ss}.spont),'--k')
        if strcmp(plot_type{ss},'User')
%             xticks(stim_label{ss}{p});
%             xticklabels(stim_name{ss})
%             xtickangle(45)
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
        subplot(2,length(SUrate{1}),ss+length(SUrate{1}))
        
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
