%% Stimulus selectivity using user call type

ss = 2;


N_list{ss} = unique([Pool{ss}.neuron_nb]);
N_list{ss} = N_list{ss}(2:end);
 
stim_list = [1,2,3,4,7,8];
%1,2 are Phees
%3,4 are Trills
%7,8 are Twitters

for c = 1:8
resp_cat{c} = [];
end

for n = 1:length(N_list{ss})
    rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list{ss}(n));
    p = rec_list{ss}(1);
    N_stim = size(rate{ss}.stim{p},1);
    stim_av = mean(rate{ss}.stim{p},2);
    pre_std = std2(rate{ss}.pre{p});
    pre_av = mean2(rate{ss}.pre{p});
    
    for st = 1:length(stim_list)
        check(st) = stim_av(stim_list(st))>pre_av +pre_std/sqrt(length(stim_list(st)));
    end
    check2 = zeros(1,3); %doesn't account for reverse
%     if sum(check)>1
%         check
%         pause
%     end
    if check(1) || check(2)    %either Phees
        check2(1) = 1;
    end
    if check(3) || check(4)    %either Trills
        check2(2) = 1;
    end
    if check(5) || check(6)    %either twitters
        check2(3) = 1;
    end
    
    
    % account for selectivity
    if check2 == [1 0 0]   %either Phees
        resp_cat{1} = [resp_cat{1}; N_list{ss}(n)];
    elseif check2 == [0 1 0]    %either Trills
        resp_cat{2} = [resp_cat{2}; N_list{ss}(n)];
    elseif check2 == [0 0 1]    %either Twitter
        resp_cat{3} = [resp_cat{3}; N_list{ss}(n)];
    elseif check2 == [1 1 0]
        resp_cat{4} = [resp_cat{4}; N_list{ss}(n)];
    elseif check2 == [1 0 1]
        resp_cat{5} = [resp_cat{5}; N_list{ss}(n)];
    elseif check2 == [0 1 1]
        resp_cat{6} = [resp_cat{6}; N_list{ss}(n)];
    elseif check2 == [1 1 1]
        resp_cat{7} = [resp_cat{7}; N_list{ss}(n)];
    else
        resp_cat{8} = [resp_cat{8}; N_list{ss}(n)];
    end
    
end
    
for c = 1:8
    cat_count(c) = length(resp_cat{c});
end
    
figure
bar(cat_count)
xticklabels([{'P'},{'Tr'},{'Tw'},{'P+Tr'},{'P+Tw'},{'Tw+Tr'},{'All'},{'None'}]);

%%
% Now use all three stimulus sets for vocalizations

N_list_all = []; 
for ss = [3 4 5] %twitter Trill Phee
    N_list{ss} = unique([Pool{ss}.neuron_nb]);
    N_list{ss} = N_list{ss}(2:end);
    N_list_all =[N_list_all unique([Pool{ss}.neuron_nb])];
end
N_list_all = unique(N_list_all);
N_list_all = N_list_all(2:end);

for c = 1:8
resp_cat2{c} = [];
end

vocal_dist = {};
vocal_dist{3} = [];
vocal_dist{4} = [];
vocal_dist{5} = [];

N_list_int = intersect(N_list{3},N_list{4});
N_list_int = intersect(N_list_int,N_list{5});


for n = 1:length(N_list_int)
    check2 = zeros(1,3);
    
    for ss = 3:5
        m = find((N_list{ss} == N_list_all(n)));
        if ~isempty(m)
            rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list{ss}(m));
            p = rec_list{ss}(1);
            N_stim = size(rate{ss}.stim{p},1);
            stim_av = mean(rate{ss}.stim{p},2);
            pre_std = std2(rate{ss}.pre{p});
            pre_av = mean2(rate{ss}.pre{p});
            
            stim_list = Pool{ss}(p).xb.stimulus_ch1(:,1);
            for st = 1:length(stim_list)
                check = stim_av(stim_list)>pre_av+pre_std/sqrt(length(stim_list));
            end
            vocal_dist{ss} = [vocal_dist{ss}; sum(check)];
            nat_list = stim_list;
%             nat_list = 1:3:20;
%             if ss == 3
%                 nat_list = [1,3,5,8,10,11,13,15,17,19];
%             end

            if sum(check(nat_list))>0
                if ss ==3 %twitter
                    check2(3) = 1;
                elseif ss ==4 %trill
                    check2(2) = 1;
                elseif ss == 5 %Phee
                    check2(1) = 1;
                end
            end
        end
    end
%     check2 = check2(3:end);
        % account for selectivity
    if check2 == [1 0 0]   %either Phees
        resp_cat2{1} = [resp_cat2{1}; N_list{ss}(m)];
    elseif check2 == [0 1 0]    %either Trills
        resp_cat2{2} = [resp_cat2{2}; N_list{ss}(m)];
    elseif check2 == [0 0 1]    %either Twitter
        resp_cat2{3} = [resp_cat2{3}; N_list{ss}(m)];
    elseif check2 == [1 1 0]
        resp_cat2{4} = [resp_cat2{4}; N_list{ss}(m)];
    elseif check2 == [1 0 1]
        resp_cat2{5} = [resp_cat2{5}; N_list{ss}(m)];
    elseif check2 == [0 1 1]
        resp_cat2{6} = [resp_cat2{6}; N_list{ss}(m)];
    elseif check2 == [1 1 1]
        resp_cat2{7} = [resp_cat2{7}; N_list{ss}(m)];
    else
        resp_cat2{8} = [resp_cat2{8}; N_list{ss}(m)];
    end
end
    
% 
% for n = 1:length(N_list_all)
%     check2 = zeros(1,3);
%     
%     for ss = 3:5
%         m = find((N_list{ss} == N_list_all(n)));
%         if ~isempty(m)
%             rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list{ss}(m));
%             p = rec_list{ss}(1);
%             N_stim = size(rate{ss}.stim{p},1);
%             stim_av = mean(rate{ss}.stim{p},2);
%             pre_std = std2(rate{ss}.pre{p});
%             pre_av = mean2(rate{ss}.pre{p});
%             
%             stim_list = Pool{ss}(p).xb.stimulus_ch1(:,1);
%             for st = 1:length(stim_list)
%                 check = stim_av(stim_list)>pre_av+pre_std/sqrt(length(stim_list));
%             end
%             vocal_dist{ss} = [vocal_dist{ss}; sum(check)];
%             nat_list = stim_list;
% %             nat_list = 1:3:20;
% %             if ss == 3
% %                 nat_list = [1,3,5,8,10,11,13,15,17,19];
% %             end
% 
%             if sum(check(nat_list))>0
%                 if ss ==3 %twitter
%                     check2(3) = 1;
%                 elseif ss ==4 %trill
%                     check2(2) = 1;
%                 elseif ss == 5 %Phee
%                     check2(1) = 1;
%                 end
%             end
%         end
%     end
% %     check2 = check2(3:end);
%         % account for selectivity
%     if check2 == [1 0 0]   %either Phees
%         resp_cat2{1} = [resp_cat2{1}; N_list{ss}(m)];
%     elseif check2 == [0 1 0]    %either Trills
%         resp_cat2{2} = [resp_cat2{2}; N_list{ss}(m)];
%     elseif check2 == [0 0 1]    %either Twitter
%         resp_cat2{3} = [resp_cat2{3}; N_list{ss}(m)];
%     elseif check2 == [1 1 0]
%         resp_cat2{4} = [resp_cat2{4}; N_list{ss}(m)];
%     elseif check2 == [1 0 1]
%         resp_cat2{5} = [resp_cat2{5}; N_list{ss}(m)];
%     elseif check2 == [0 1 1]
%         resp_cat2{6} = [resp_cat2{6}; N_list{ss}(m)];
%     elseif check2 == [1 1 1]
%         resp_cat2{7} = [resp_cat2{7}; N_list{ss}(m)];
%     else
%         resp_cat2{8} = [resp_cat2{8}; N_list{ss}(m)];
%     end
% end
%     
for c = 1:8
    cat_count(c) = length(resp_cat2{c});
end

% cat_count(8) = cat_count(8) + length(N_list_all)-sum(cat_count);

sum(cat_count);
figure
bar(cat_count)
xticklabels([{'P'},{'Tr'},{'Tw'},{'P+Tr'},{'P+Tw'},{'Tw+Tr'},{'All'},{'None'}]);

figure

subplot(1,3,1)
histogram(vocal_dist{3},[0:1:20])
title('twitters')
subplot(1,3,2)
histogram(vocal_dist{4},[0:1:20])
title('trills')
subplot(1,3,3)
histogram(vocal_dist{5},[0:1:20])
title('phees')




%%  
    