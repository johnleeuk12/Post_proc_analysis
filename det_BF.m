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
                switch db
                    case 60
                        s2 = 1;
                    case 40
                        s2 = 2;
                    case 20
                        s2 = 3;
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
                StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
                %             stim_label{ss}{p} = Pool{ss}(p).xb.stimulus_ch1(:,8);
            end
        end
    end
end

        
    
