function Quality_control(Pool, rate, raster)

%% comparing PT first and PT last

N_pool = size(Pool,2); %ss will be 1 in this case
ss = 1;
N_list{ss} = unique([Pool{ss}.neuron_nb]);
N_list{ss} = N_list{ss}(2:end);


h = 5;
D_param = zeros(1,length(N_list{ss}));
D_param2 = zeros(1,length(N_list{ss}));
tic
for n = 1:length(N_list{ss}) % n is each neuron
    rec_list = {};
    check_counter = 1;
%     SUrate{n} = {};
    
    rec_list{ss} = find([Pool{ss}.neuron_nb] ==N_list{ss}(n));
    if isempty(rec_list{ss})
        check_counter = 0;
    elseif length(rec_list{ss})<2
        check_counter = 0;
    end
    
    if check_counter ~= 0
        p(1) = rec_list{ss}(1);
        p(2) = rec_list{ss}(end);
        av_PSTH = {};
        smooth_PSTH ={};
        smooth_PSTH{1} = [];
        smooth_PSTH{2} = []; 
        av_PSTH{1} = [];
        av_PSTH{2} = [];
        for m = 1:2
            
            N_stim = Pool{ss}(p(m)).xb.nstim;
            
%             figure(m)
            D = length(rate{ss}.PSTH{p(m),1});
            for st = 1: N_stim
                xs = 1:D;
                av_PSTH{m} = [av_PSTH{m}; mean(rate{ss}.PSTH{p(m),st},1)];
                for i = 1:D
                    ys(i) = gaussian_kern_reg(xs(i),xs,mean(rate{ss}.PSTH{p(m),st},1),h);
                end
                smooth_PSTH{m} = [smooth_PSTH{m}; ys];
%                 plot(ys)
%                 hold on
%                 drawnow
            end
%             hold off
            
            spont_rate(m) = mean(mean(rate{ss}.pre{p(m)}));
            smooth_PSTH{m} = smooth_PSTH{m}-spont_rate(m);
            
        end  
        % Cross validation of euclidean distance
%         D_temp = zeros(1,100);
%         for t = 1:100
%             ind1 = randsample(1:N_stim,round(N_stim/2));
%             ind2 = setdiff(1:N_stim,ind1);
%             temp1 = zeros(size(smooth_PSTH{1}));
%             temp2 = zeros(size(smooth_PSTH{1}));
%             
%             temp1(ind1,:) = smooth_PSTH{1}(ind1,:);
%             temp1(ind2,:) = smooth_PSTH{2}(ind2,:);
%             temp2(ind1,:) = smooth_PSTH{2}(ind1,:);
%             temp2(ind2,:) = smooth_PSTH{1}(ind2,:);
%             
%             max_FR = max(max(max(abs(temp1))), max(max(abs(temp2))));
%             
%             temp1 = temp1/max_FR;
%             temp2 = temp2/max_FR;
%             
%             
%             D_temp(t) = norm(temp2-temp1);
%         end
%         
%         D_rand = mean(D_temp);
        
        
        max_FR = max(max(max(abs(smooth_PSTH{2}))), max(max(abs(smooth_PSTH{1}))));
        for m =1:2
            smooth_PSTH{m} = smooth_PSTH{m}/max_FR;
        end
        D_param(n) = norm(smooth_PSTH{1}- smooth_PSTH{2}); %-D_rand;
        

        
        if D_param(n)<7
            figure(n)
            sgtitle(['D = ', num2str(D_param(n))])
            subplot(1,2,1)
            imagesc(smooth_PSTH{1})
            %         s = surf(smooth_PSTH{1});
            %         s.EdgeColor = 'none';
            subplot(1,2,2)
            imagesc(smooth_PSTH{2})
            
            %         s2 = surf(smooth_PSTH{2});
            %         s2.EdgeColor = 'none';
            drawnow
        end
    end
    if D_param(n) ==0
        D_param(n) = nan;
    end
    toc
end


for n = 1:127
        if D_param(n) ==0
        D_param(n) = nan;
        end
end

% list = find(D_param>7);

