function DD = format4DataHigh(SUrate,vocal_type)

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
vocal_type = 31;

switch vocal_type
    case 11
        vc_list1 = round([7.16:0.12:9.56]*100)/100;
        stim_dur = 1180;
    case 12
        vc_list1 = round([4.76:0.12:7.16]*100)/100;
        stim_dur = 1180;

    case 21
        vc_list1 = round([6.64:0.41:10.74]*100)/100;
        vc_fmr = 29.96; % change here to compare
        stim_dur = 406;
    case 22
        vc_list1 = round([2.95:0.41:6.64]*100)/100;
        vc_fmr = 29.96; % change here to compare
        stim_dur = 406;
    case 31
        vc_list1 = round([9.72:0.735:17.07]*100)/100;
        vc_ipi = 139.8;
        stim_dur = 1161;
    case 32
        vc_list1 = round([2.37:0.735:9.72]*100)/100;
        vc_ipi = 0.1398;
end





% list of vocalization parameters to chose from. Change parameters here
% trill CF  %
% vc_list1 = round([6.64:0.41:10.74]*100)/100;
% vc_list1 = round([2.95:0.41:6.64]*100)/100;
% vc_list1 = round([2.95:0.41:6.64]*100)/100;
% vc_list2 = [25,32];
% vc_list3 = [900, 1000];
% 




NN = length(SUrate); 

% find list of neurons that responded to all stim on list. \
good_list = [];
good_list_pid = [];
good_list_nid = [];
for n = 1:NN
    try
        if intersect(round((unique(SUrate{n}{1}.stim(:,1)).')*1e2),round(vc_list1*1e2)) == round(vc_list1*1e2)
            switch vocal_type
                case {01,02, 11,12,31, 32}
                    good_list = [good_list,n];
                    good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                    good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                case {21, 22}
                    if length(unique(SUrate{n}{1}.stim(:,3))) ==1 % fm depth has 1 value
                        good_list = [good_list,n];
                        good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                        good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                    end
                case {23,24}
                    if length(unique(SUrate{n}{1}.stim(:,2))) ==1 % fm rate  has 1 value
                        good_list = [good_list,n];
                        good_list_pid = [good_list_pid, SUrate{n}{1}.pid];
                        good_list_nid = [good_list_nid, SUrate{n}{1}.nid];
                    end
                
            end
            
        end
    catch
    end
    
    
end

%%

DD= {};
nreps = 8;

pre_stim = 300;
post_stim = 500;
cmap1 = parula(length(vc_list1));
cmap2 = gray(length(vc_list1));


for p = 1:length(vc_list1)
    vc_cf = vc_list1(p);
    
    for r = 1:nreps

        trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
        DD(r+(p-1)*nreps).data = zeros(length(good_list),trial_dur);
        DD(r+(p-1)*nreps).condition = num2str(vc_cf);
        DD(r+(p-1)*nreps).epochStarts = [1,pre_stim,pre_stim+stim_dur];
        DD(r+(p-1)*nreps).epochColors = [0.7, 0.7, 0.7;cmap1(p,:);cmap2(p,:)];
        for n = 1:length(good_list)
            switch vocal_type
                case {01, 02}
                    vc_ind = find(round(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,8)) == round(vc_cf*1e3)); 
                case {11, 12}
                    vc_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == vc_cf);
                case {21, 22}
                    vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_fmr);
                case {31, 32}
                    vc_ind = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf & SUrate{good_list(n)}{1}.stim(:,2) == vc_ipi);
            end
            
            DD(r+(p-1)*nreps).data(n,:) = round(SUrate{good_list(n)}{1}.PSTH{vc_ind}(r,1:trial_dur)*1e-3);
        end
    end
end



%%






















