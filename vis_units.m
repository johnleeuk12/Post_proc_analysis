function vis_units(list)

% function to visualize individual neurons 
list = flip(new_SU_list);
N = length(list);

%% 


%%
for n = 530:540
    figure(n)
    
    for t = 1:2
        raster = {};
        raster.spikes = [];
        raster.stim = [];
        raster.reps = [];
        
        %     PSTH = [];
        for st = 1:stim_num
            for r = 1:nreps
                
                %             spikes = find(SUrate{list(n)}{1}.PSTH{st}(r,:)>0);
                
                if t ==1
                    spikes = find(D_phee((st-1)*nreps+r).data(n,:)>0);
                else
                    spikes = find(D_PT((st-1)*nreps+r).data(n,:)>0);
                end
                
                raster.spikes =[raster.spikes spikes];
                raster.reps = [raster.reps ones(1,length(spikes))*r];
                raster.stim = [raster.stim ones(1,length(spikes))*st];
            end
        end
        raster.spikes = raster.spikes-PreStim;
        
        
        subplot(1,2,t)
        rectangle('Position',[0 0,stim_dur nreps*stim_num],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
        hold on
        
        plot(raster.spikes,nreps*(raster.stim-1)+raster.reps,'k.','MarkerSize',15);
        yticks([1:nreps*2:nreps*stim_num]+5)
        stim_name = {};
        for st = 1:stim_num
            stim_name{st} = num2str(stim_set(st));
        end
        yticklabels(stim_name(1:2:end));
        xlabel('time (s)')
        axis([-PreStim max(stim_dur) + PostStim 0 nreps*stim_num+1])
        set(gcf, 'Position', [800 800 1000 400]); 
    end
    if na_id_both(2,n) == 1
        sgtitle(['M60Fu' num2str(na_id_both(1,n))])
    elseif na_id_both(2,n) == 2
        sgtitle(['M160Eu' num2str(na_id_both(1,n))])
    end
    pause(0.1)
    
end