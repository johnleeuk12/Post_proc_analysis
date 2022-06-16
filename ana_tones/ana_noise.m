function out =ana_noise(H_nb_list,animal_name,B,N_list2,t_end)
%% analysis comparing PT and noise 
%positive values for percentage indicate an increase for noise. 

% 
% H_nb_list = 1;
% animal_name = 'M56E';
noise_ind = 2;

out = {};
out.rawdata = {};
out.data = [];
out.info = [];


savedir = fullfile('E:\DATA',filesep,animal_name,filesep,'ana_tones\data');
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
N_list = [neurons_loc_tag.neuron_nb];



tic
nn =1;
for H_nb = H_nb_list
    hole_ind = find([neurons_loc_tag.hole_nb] == H_nb);
    tracks = unique([neurons_loc_tag(hole_ind).track_nb]);
    
    % figure
    for t = 1:length(tracks)
        fprintf(['processing H' num2str(H_nb) 'T' num2str(tracks(t)) '     time : %6.2f sec \n'],toc')
        
        sub_n_list = N_list(find([neurons_loc_tag.track_nb] ==tracks(t) &[neurons_loc_tag.hole_nb] == H_nb ));
        for n =1:length(sub_n_list)
            nid = sub_n_list(n);
            PTname= 'PTn0000';
            PTname = [PTname(1:end-length(num2str(nid))) num2str(nid) '.mat'];
            
            try
                load(PTname);
                try
                    ana_PT = filename; %currently loaded file is filename. change afterwards
                catch
                end
                
                
                st_mat = zeros(1,12);
                for st_ind = 1:12
                    try
                        st_mat(st_ind) = ~isempty(ana_PT{1,st_ind});
                    catch
                    end
                end
                st_mat = reshape(st_mat,4,3);
                db_ind = find(mean(st_mat(:,1:noise_ind),2) == 1);
                st_ind2 = 1;
                
                if size(ana_PT,1) ==2 && isempty(db_ind)
                    st_mat = zeros(1,12);
                    for st_ind = 1:12
                        try
                            st_mat(st_ind) = ~isempty(ana_PT{2,st_ind});
                        catch
                        end
                    end
                    st_mat = reshape(st_mat,4,3);
                    
                    db_ind = find(mean(st_mat(1:noise_ind),2) == 1);
                    st_ind2 = 2;
                    
                end
                
                
                % change here if want to use bw = 0.5
                %                 if length(ana_PT) > 8
                %                     if ~isempty(ana_PT{st_ind2,db_ind+noise_ind*4})
                %                     noise_ind = noise_ind+1;
                %                     end
                %                 end
                
                
                %extracting response data
                %                 figure
                
                if length(db_ind) > 1
                    db_ind = B(find(N_list2 == nid),1);
                end
                
                X = {};
                for nz = 1:noise_ind
                    
                    % case where only 1 dB has PT and noise.
                    if length(db_ind) == 1
                        stim_freq = ana_PT{st_ind2,db_ind+(nz-1)*4}.xb.stimulus_ch1(:,8);
                        
                        X{nz} = [];
                        for st = 1:length(stim_freq)
                            X{nz} = [X{nz}; mean(ana_PT{st_ind2,db_ind+(nz-1)*4}.PSTH{st},1)]; %-mean(ana_PT{db_ind}.spont)]
                        end
                    else
                    end
                    out.rawdata{nn,nz} = X{nz}; % saving raw data
                    
                    % finding peak and BF
                    X2 = imgaussfilt(X{nz},[3,20],'Padding',0); %,'Padding','circular');
                    X2 = X2-mean2(X2(:,1:150)); % removing spont rate
%                     subplot(1,2,nz)
                    %                     imagesc(X2)
                    % change t_end if we're defining search window
%                     t_end = 350;
                    
                    [BF_ind, pos, ~]  = find(X2 == max(max(X2(:,200:t_end))));
                    
                    BF_ind_keep = [];
                    for bff = 1:length(BF_ind)
                        if X2(BF_ind(bff),pos(bff)) > std2(X2(:,1:150))*4 % change here for BF detection threshold
                            BF_ind_keep = [BF_ind_keep,bff];
                        end
                    end
                    BF_ind = BF_ind(BF_ind_keep);
                    pos = pos(BF_ind_keep);
                    
                    % for now, remove responses with multiple peaks of the
                    % same FR ( unlikely)
                    if length(BF_ind)>1
                        BF_ind = nan;
                    end
                    if pos < 203 %if peak is found during the first 2ms of onset
                        BF_ind = [];
                    end
                    
                    if isempty(BF_ind) % 0 means that data exists but no peak was found
                        BF_ind = 0;
                    end
                    

                    if BF_ind ~=0
                        out.data(nn,nz+3) = mean(X2(BF_ind, pos-25:pos+25)); % FR calculated within a 50ms bin centered at peak position
                    else
                        out.data(nn,nz+3) = 0;
                    end
                    
                    if  mean(X2(BF_ind, pos-25:pos+25))<1
                         out.data(nn,nz+3) = 0;
                         BF_ind = 0;
                    end

                    out.data(nn,1) = nid;
                    out.data(nn,nz+1) = BF_ind;
                end
                
                out.info(nn,1) = nid;
                out.info(nn,2) = H_nb;
                out.info(nn,3) = tracks(t);
                


                if out.data(nn,2) ==0
                    out.data(nn,6) = (out.data(nn,5)-out.data(nn,4))/max(out.data(nn,4:5));
                elseif out.data(nn,3) ==0
                    out.data(nn,6) = (out.data(nn,5)-out.data(nn,4))/max(out.data(nn,4:5));
                elseif abs(out.data(nn,2)-out.data(nn,3)) <6
                    out.data(nn,6) = (out.data(nn,5)-out.data(nn,4))/max(out.data(nn,4:5));
                else
                    out.data(nn,6) = nan;
                       
                end
                if out.data(nn,2) == 0 && out.data(nn,3) ==0
                else
                    nn = nn +1;
                end
                
                % calculating change in FR for units with similar responses
                % between PT and noise

            catch
                
            end

        end
    end
end
    

% list = find(~isnan(out.data(:,6)));
% figure
% scatter(out.data(list,4),out.data(list,5))
% % hold on
% 
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
