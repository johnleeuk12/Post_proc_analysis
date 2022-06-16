function [DD, good_list_nid,good_list,good_list_aid] = format4DataHigh(SUrate,vocal_type,vocal_param2)

%{
D = 1xN struct array with fields:
    data: (number of neurons x number of milliseconds) array
    N is the number of trials..
%}



%% Initialize, pooling neurons, checking number of trials
%{
PTPhee High         01
PTPhee Low          02

PheeCF High         11
PheeCF Low          12

TrillCF High        21
TrillCF Low         22
             
TwitterCF High      31
TwitterCF Low       32
          
%}
% vocal_type = 32;


NN = length(SUrate);
animal_list = {'M60F','M160E','M56E'};


switch vocal_type
    
    case 0
        vc_list1 = round([4.76:0.12:9.56]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,8)*1e-3;
        end
        
        
    case 01
        vc_list1 = round([6.68:0.12:9.56]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,8)*1e-3;
        end
    case 02
        vc_list1 = round([4.76:0.12:7.64]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,8)*1e-3;
        end
    case 10
        vc_list1 = round([4.76:0.12:9.56]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,10);
        end
        
    case 11
        vc_list1 = round([6.68:0.12:9.56]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,10);
        end
    case 12
        vc_list1 = round([4.76:0.12:7.64]*100)/100;
        stim_dur = 1180;
        nreps = 10;
        
        for n = 1:NN
            SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,10);
        end
    case 21
        vc_list1 = round([5.82:0.41:10.74]*100)/100;
        vc_fmr_list = round([29.96:2.69:56.86]*100)/100;
        vc_fmr_ind = vocal_param2;
        vc_fmr = vc_fmr_list(vc_fmr_ind); %1 to 11
        vc_list2 = vc_fmr_list(vocal_param2); %: vocal_param2+2);
        
        %                         vc_fmr = 29.96; % change here to compare
        %         vc_fmr = 32.65; %0.5SD
        %         vc_fmr = 35.34; %1SD
        %         vc_fmr = 38.03; % 1.5SD
        %         vc_fmr = 40.72; %2SD
        %         vc_fmr = 43.41; % 2.5SD
        %         vc_fmr = 51.48; %4SD
        stim_dur = 406;
        nreps = 10;
        
    case 22
        vc_list1 = round([2.54:0.41:7.46]*100)/100;
        vc_fmr_list = round([29.96:2.69:56.86]*100)/100;
        vc_fmr_ind = vocal_param2;
        vc_fmr = vc_fmr_list(vc_fmr_ind); %1 to 11
        vc_list2 = vc_fmr_list(vocal_param2);% : vocal_param2+2);
        %                 vc_fmr = 29.96; % change here to compare
        %         vc_fmr = 32.65; %0.5SD
        %         vc_fmr = 35.34; %1SD
        %         vc_fmr = 38.03; % 1.5SD
        %         vc_fmr = 40.72; %2SD
        %         vc_fmr = 43.41; % 2.5SD
        %         vc_fmr = 51.48; %4SD
        stim_dur = 406;
        nreps = 10;
        
    case 31
        vc_list1 = round([8.25:0.735:17.07]*100)/100;
        %         vc_ipi = 0.1580;
        %         stim_dur = 1306;
        %
        %                 vc_ipi = 0.2309;
        %                 stim_dur = 1890;
        vc_ipi_list = [0.1398, 0.1489, 0.1580, 0.1671, 0.1763, 0.1854, 0.1945, 0.2036, 0.2127, 0.2218, 0.2309];
%         vc_list2 = vc_ipi_list(vocal_param2:vocal_param2+2);
        vc_list2 = vc_ipi_list(vocal_param2);
        stim_dur_list = [1161, 1233:73:1890];

        stim_dur = stim_dur_list(vocal_param2);
        nreps = 8;
        
    case 32
        vc_list1 = round([2.37:0.735:11.19]*100)/100;
        vc_ipi_list = [0.1398, 0.1489, 0.1580, 0.1671, 0.1763, 0.1854, 0.1945, 0.2036, 0.2127, 0.2218, 0.2309];
%         vc_list2 = vc_ipi_list(vocal_param2:vocal_param2+2);
        vc_list2 = vc_ipi_list(vocal_param2);
        stim_dur_list = [1161, 1233:73:1890];

        stim_dur = stim_dur_list(vocal_param2);
        nreps = 8;
        
end












% find list of neurons that responded to all stim on list. \
good_list = [];
good_list_pid = [];
good_list_nid = [];
good_list_aid = [];
for n = 1:NN
    try
        if intersect(round((unique(SUrate{n}{1}.stim(:,1)).')*1e2),round(vc_list1*1e2)) == round(vc_list1*1e2)
            switch vocal_type
                case {0 01,02,10,11,12}
                    good_list = [good_list,n];
                    good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                    good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                    for a = 1:length(animal_list)
                        if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
                            good_list_aid = [good_list_aid a];
                        end
                    end
                case {21, 22}
                    if length(unique(SUrate{n}{1}.stim(:,3))) ==1 % fm depth has 1 value
                        %                         if intersect(round((unique(SUrate{n}{1}.stim(:,2)).')*1e2),round(vc_fmr_list*1e2)) == round(vc_fmr_list*1e2)
                        if size(SUrate{n}{1}.raw,2) >= nreps
                            good_list = [good_list,n];
                            good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                            good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                            for a = 1:length(animal_list)
                                if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
                                    good_list_aid = [good_list_aid a];
                                end
                            end
%                         else
%                             good_list = [good_list,n];
%                             good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
%                             good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
%                             for a = 1:length(animal_list)
%                                 if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
%                                     good_list_aid = [good_list_aid a];
%                                 end
%                             end
%                             for nst = 1:length(SUrate{n}{1}.PSTH)
%                                 SUrate{n}{1}.PSTH{nst} = [SUrate{n}{1}.PSTH{nst} ;SUrate{n}{1}.PSTH{nst}];
%                             end
                        end
                        
                        %                         end
                    end
                case {23,24}
                    if length(unique(SUrate{n}{1}.stim(:,2))) ==1 % fm rate  has 1 value
                        good_list = [good_list,n];
                        good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                        good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                        for a = 1:length(animal_list)
                            if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
                                good_list_aid = [good_list_aid a];
                            end
                        end
                    end
                    % case for m60F
                    %                 case {31,32}
                    % %                     T1 = SUrate{n}{1}.spont(1:66);
                    % %                     T2 = SUrate{n}{1}.spont(78:end);
                    % %                     if ttest2(T1,T2) == 0 % problem with spont rate in twitters
                    %                         good_list = [good_list,n];
                    %                         good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                    %                         good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                    % %                     end
                case{31,32}
                    if size(SUrate{n}{1}.PSTH{1},1) > 7
                        if length(SUrate{n}{1}.spont) >1
                            T1 = SUrate{n}{1}.spont(1:66);
                            T2 = SUrate{n}{1}.spont(78:end);
                            if ttest2(T1,T2) == 0 % problem with spont rate in twitters
                                good_list = [good_list,n];
                                good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                                good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                                for a = 1:length(animal_list)
                                    if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
                                        good_list_aid = [good_list_aid a];
                                    end
                                end
                            end
                        else
                            good_list = [good_list,n];
                            good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                            good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                            for a = 1:length(animal_list)
                                if strcmp(SUrate{n}{1}.xb.animal,animal_list(a))
                                    good_list_aid = [good_list_aid a];
                                end
                            end
                        end
                    end
            end
            
        end
    catch
    end
    
    
end

%%

DD= {};
% nreps = 8;

pre_stim = 300;
post_stim = 500;
cmap1 = parula(length(vc_list1));
cmap2 = gray(length(vc_list1));

switch vocal_type
    case {0, 01, 02, 10, 11, 12}
        for p = 1:length(vc_list1)
            vc_cf = vc_list1(p);
            
            for r = 1:nreps
                
                trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
                DD(r+(p-1)*nreps).data = zeros(length(good_list),trial_dur);
                DD(r+(p-1)*nreps).condition = num2str(vc_cf,'%05.2f');
                DD(r+(p-1)*nreps).epochStarts = [1,pre_stim,pre_stim+stim_dur];
                DD(r+(p-1)*nreps).epochColors = [0.7, 0.7, 0.7;cmap1(p,:);cmap2(p,:)];
                for n = 1:length(good_list)
                    switch vocal_type
                        case {0, 01, 02}
                            vc_ind = find(round(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,8)) == round(vc_cf*1e3));
                        case {10,11, 12}
                            vc_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == vc_cf);
                        case {21, 22}
                            vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_fmr_list(vc_fmr_ind));
                            vc_ind2 = vc_ind+1;
                            %                     vc_ind2 = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_fmr_list(vc_fmr_ind+1));
                            if isempty(vc_ind)
                                [~,vc_ind]=min(abs(SUrate{good_list(n)}{1}.stim(:,2)-vc_fmr));
                                vc_ind = vc_ind+(p-1)*length(unique(SUrate{good_list(n)}{1}.stim(:,2)));
                                vc_ind2 = vc_ind +1;
                            end
%                         case {31, 32}
%                             SUrate{good_list(n)}{1}.stim(:,1) = round(SUrate{good_list(n)}{1}.stim(:,1)*1e2)/1e2;
%                             vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_ipi);
                    end
                    
                    if vocal_type == 21 || vocal_type == 22 % attempting to group 0 ~1SD together
                        PSTH_group = (SUrate{good_list(n)}{1}.PSTH{vc_ind}(r,1:trial_dur) + SUrate{good_list(n)}{1}.PSTH{vc_ind2}(r,1:trial_dur))/2;
                        DD(r+(p-1)*nreps).data(n,:) = round(PSTH_group*1e-3);
                    else
                        DD(r+(p-1)*nreps).data(n,:) = round(SUrate{good_list(n)}{1}.PSTH{vc_ind}(r,1:trial_dur)*1e-3);
                    end
                end
            end
        end
        
        
    case{21,22}
        for p = 1:length(vc_list1)
            for v = 1:length(vc_list2)
                vc_cf = vc_list1(p);
                vc_fmr = vc_list2(v);
                for r = 1:nreps
                    trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).data = zeros(length(good_list),trial_dur);
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).condition = num2str(vc_cf,'%05.2f');
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).epochStarts = [1,pre_stim,pre_stim+stim_dur];
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).epochColors = [0.7, 0.7, 0.7;cmap1(p,:);cmap2(p,:)];
                    for n = 1:length(good_list)
                        vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_fmr);
                        if isempty(vc_ind)
                            [~,vc_ind]=min(abs(SUrate{good_list(n)}{1}.stim(:,2)-vc_fmr));
                            vc_ind = vc_ind+(p-1)*length(unique(SUrate{good_list(n)}{1}.stim(:,2)));
                        end
                        PSTH_group = (SUrate{good_list(n)}{1}.PSTH{vc_ind}(r,1:trial_dur) + SUrate{good_list(n)}{1}.PSTH{vc_ind+1}(r,1:trial_dur))/2;
                        DD(r+(p-1)*nreps*length(vc_list2)).data(n,:) = round(PSTH_group*1e-3);
                        
                        
                    end
                end
            end
        end
        
    case {31, 32}
        for p = 1:length(vc_list1)
            for v = 1:length(vc_list2)
                vc_cf = vc_list1(p);
                vc_ipi = vc_list2(v);
                
                for r = 1:nreps
                    
                    trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).data = zeros(length(good_list),trial_dur);
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).condition = num2str(vc_cf,'%05.2f');
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).epochStarts = [1,pre_stim,pre_stim+stim_dur];
                    DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).epochColors = [0.7, 0.7, 0.7;cmap1(p,:);cmap2(p,:)];
                    for n = 1:length(good_list)
                        
                        SUrate{good_list(n)}{1}.stim(:,1) = round(SUrate{good_list(n)}{1}.stim(:,1)*1e2)/1e2;
                        vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_ipi);
                        DD(r+(p-1)*nreps*length(vc_list2) + (v-1)*nreps).data(n,:) = round(SUrate{good_list(n)}{1}.PSTH{vc_ind}(r,1:trial_dur)*1e-3);
                    end
                end
                
                
                
            end
        end
        
        
end




%%






















