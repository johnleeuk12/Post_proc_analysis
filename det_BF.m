function det_BF(Pool, rate, raster)


%%

N_list = unique([Pool.neuron_nb]);
N_list = N_list(2:end);
SUrate = {};

for n = 1:length(N_list)
    rec_list = find([Pool.neuron_nb] == N_list(n));
    check_counter = 1;
    SUrate{n} = {};
    
    for r = 1: length(rec_list)
        p= rec_list(r);
        a_code = Pool(p).xb.analysis_code;
        db = Pool(p).xb.stimulus_ch1(1,3);
        lens = Pool(p).xb.stimulus_ch1(1,5);
        if db ~= 80
            if lens ==100
                if a_code == 1
                    switch db
                        case 60
                            s2 = 1;
                        case 40
                            s2 = 2;
                        case 20
                            s2 = 3;
                    end
                elseif a_code == 62
                    switch db
                        case 60
                            s2 = 4;
                        case 40
                            s2 = 5;
                        case 20
                            s2 = 6;
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
                SUrate{n}{s1,s2}.error = std(rate.stim{p},1,2)/sqrt(nreps);
                SUrate{n}{s1,s2}.spont = mean(rate.pre{p},2);
                SUrate{n}{s1,s2}.raw = rate.stim{p};
                SUrate{n}{s1,s2}.PSTH = {};
                SUrate{n}{s1,s2}.nid = Pool(p).neuron_nb;
                SUrate{n}{s1,s2}.pid = p;
                SUrate{n}{s1,s2}.xb = Pool(p).xb;
                SUrate{n}{s1,s2}.xb.data = [];
                
                for te = 1:size(rate.PSTH,2)
                    SUrate{n}{s1,s2}.PSTH{te}= rate.PSTH{p,te};
                end
%                 StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
                %             stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
            end
        end
    end
end

%% plotting figures

for n = 1:length(N_list)
    fi = figure(n);
    set(fi, 'Position', [400 100 1700 800]);
    for s1 =1:2
        for s2 =1:6
            try if ~isempty(SUrate{n}{s1,s2})
                    subplot(2,6,(s1-1)*12+1+(s2-1)*2+1)
                    hold off
                    stim_label = SUrate{n}{s1,s2}.xb.stimulus_ch1(:,8);
                    stim_ticks = {};
                    
                    for stim = 1:length(stim_label)
                        stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                    end
                    
                    p  =SUrate{n}{s1,s2}.pid;
                    %tuning curve
                    errorbar(SUrate{n}{s1,s2}.mean,stim_label,SUrate{n}{s1,s2}.error,'horizontal','LineWidth',2);
                    hold on
                    xline(mean(SUrate{n}{s1,s2}.spont),'--k');
                    set(gca, 'YScale', 'log')
                    axis([0 inf min(stim_label) max(stim_label)])
                    % rasterplot
                    subplot(2,6,(s1-1)*12+1+(s2-1)*2)
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

                    
                    
                    

                    yticklabels(stim_ticks(1:2:end))
                    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                    hold off
                    title('rasterplot')

                    drawnow
                    pause(0.1)
                    
                end
            catch
            end
        end
    end
    
end



















