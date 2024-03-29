function SUrate = ana_2Dstim(Pool, rate, raster, figure_on, VT_name)
%%
%{
Pooling data for 2D stim or stimuli sets that run across multiple
recording segments. 
%}


%%



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
switch VT_name
    case {'trill1', 'trill2'}
        V_ind = 10;
        st_range.cf = 2.54:0.41:10.74; %range for cf
        st_range.fmr = 3.06:2.69:56.86; %range for fm rate
        st_range.fmd = 50:90:1400;
    case 'twitter'
        V_ind = 9;
end
%%


for n = 1:length(N_list_int)

   if  length(unique(Pool{1}(n).xb.stimulus_ch1(:,V_ind)))==1% || ...
            %length(unique(Pool{1}(n).xb.stimulus_ch1(:,22)))==1
        N_list_int(n) = [];
   end
end

N_list_int = N_list_int(2:end);


N_list_int = N_list_int(2:end);
rec_list = {};
tic;
for n = 1:length(N_list_int) % n is each neuron
    check_counter = 1;
    SUrate{n} = {};
    
    
    for ss = 1:N_pool %ss is each stim for same neuron
        rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list_int(n));
        if isempty(rec_list{ss})
            check_counter = 0;
            %         elseif length(rec_list{ss})<2 %This is for Puretones
            %             check_counter = 0;
        end
        
        
        %     ss = 1;
        if check_counter == 1
            SUrate{n}{ss}.mean = [];
            SUrate{n}{ss}.error = [];
            SUrate{n}{ss}.spont = [];
            SUrate{n}{ss}.raw = [];
            SUrate{n}{ss}.stim = [];
            SUrate{n}{ss}.pid = [];
            SUrate{n}{ss}.raster = {};
            SUrate{n}{ss}.raster.stim = [];
            SUrate{n}{ss}.raster.rep = [];
            SUrate{n}{ss}.raster.spikes = [];
            
            ana_code = Pool{ss}(1).xb.analysis_code;
            
            if ana_code < 2000
                rec_list{ss} = rec_list{ss}(1); % or end
            end
            
            for r = 1:length(rec_list{ss})
                p = rec_list{ss}(r);
                nreps = size(rate{ss}.stim{p},2);
                SUrate{n}{ss}.mean = [SUrate{n}{ss}.mean ; mean(rate{ss}.stim{p},2)];
                SUrate{n}{ss}.error = [SUrate{n}{ss}.error; std(rate{ss}.stim{p},1,2)/sqrt(nreps)];
                SUrate{n}{ss}.spont = [ SUrate{n}{ss}.spont; mean(rate{ss}.pre{p},2) ];
                SUrate{n}{ss}.raw = [ SUrate{n}{ss}.raw; rate{ss}.stim{p}] ;
                SUrate{n}{ss}.nid = Pool{ss}(p).neuron_nb;
                SUrate{n}{ss}.pid = [SUrate{n}{ss}.pid, p];
                SUrate{n}{ss}.xb = Pool{ss}(p).xb;
                SUrate{n}{ss}.xb.data = [];
                

                if r ==1
                    SUrate{n}{ss}.PSTH = {};
                    for te = 1:size(rate{ss}.PSTH,2)
                        if ~isempty(rate{ss}.PSTH{p,te})
                            SUrate{n}{ss}.PSTH{te}= rate{ss}.PSTH{p,te};
                        end
                    end
                    SUrate{n}{ss}.raster.stim = raster{1}.stim{p};
                    SUrate{n}{ss}.raster.rep = raster{1}.rep{p};
                    SUrate{n}{ss}.raster.spikes = raster{1}.spikes{p};
                    
                else
                    te_ind = length(SUrate{n}{ss}.PSTH);
                    for te = 1:size(rate{ss}.PSTH,2)
                        if ~isempty(rate{ss}.PSTH{p,te})
                            SUrate{n}{ss}.PSTH{te+te_ind}= rate{ss}.PSTH{p,te};
                        end
                    end
                    SUrate{n}{ss}.raster.stim = [SUrate{n}{ss}.raster.stim, raster{1}.stim{p}+Pool{1}(rec_list{ss}(1)).xb.nstim];
                    SUrate{n}{ss}.raster.rep = [SUrate{n}{ss}.raster.rep, raster{1}.rep{p}];
                    SUrate{n}{ss}.raster.spikes = [SUrate{n}{ss}.raster.spikes, raster{1}.spikes{p}];
                end
                
                switch ana_code
                    case {1,62}
                        stid_1 = 8;
                        SUrate{n}{ss}.stim_tags = {{'cf'}};
                        SUrate{n}{ss}.stim = [SUrate{n}{ss}.stim; SUrate{n}{ss}.xb.stimulus_ch1(:,stid_1)];
                    case 2335 % trills
                        stid_1 = 10;
                        stid_2 = 22;
                        stid_3 = 24;
                        SUrate{n}{ss}.stim_tags = {{'cf'},{'fmrate'},{'fmdepth'}};
                        SUrate{n}{ss}.stim = [SUrate{n}{ss}.stim; SUrate{n}{ss}.xb.stimulus_ch1(:,[stid_1, stid_2, stid_3])];
                    case 2367 % Twitter
                        stid_1 = 9;
                        stid_2 = 12;
                        stid_3 = 5;
                        SUrate{n}{ss}.stim_tags = {{'cf'},{'IPI'},{'Stim_dur'}};
                        SUrate{n}{ss}.stim = [SUrate{n}{ss}.stim; SUrate{n}{ss}.xb.stimulus_ch1(:,[stid_1, stid_2, stid_3])];
                end
            end
            
        end
    end
end
toc
% x1 = Pool{1, 1}(72).xb.stimulus_ch1;
% x2 = Pool{1, 1}(105).xb.stimulus_ch1;
% 
% cf1 = unique(x1(:,10));
% cf2 = unique(x2(:,10));
% 








