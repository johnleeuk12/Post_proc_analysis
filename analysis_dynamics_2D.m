%% PSTH
N = length(SUrate);

% list of vocalization parameters to chose from. Change parameters here
% trill CF  %
% vc_list1 = round([6.64:0.41:10.74]*100)/100;
% vc_list1 = round([2.95:0.41:6.64]*100)/100;
% vc_list1 = round([2.95:0.41:6.64]*100)/100;
% vc_list2 = [25,32];
% vc_list3 = [900, 1000];
% 
% set2 = [10.33,10.74];
% set1 = [6.64,7.05];
% test_set =  round([7.46:0.41:9.92]*100)/100;

% Phee CF 
vc_list1 = round([4.76:0.12:7.16]*100)/100;
set2 = [4.76, 4.88, 5];
set1 = [6.92,7.04, 7.16];
test_set =  round([5.12:0.12:6.8]*100)/100;



% find list of neurons that responded to all stim on list. \
good_list = [];
for n = 1:N
    SUrate{n}{1}.stim = SUrate{n}{1}.xb.stimulus_ch1(:,10);
    try
        if intersect(unique(SUrate{n}{1}.stim(:,1)).',vc_list1) == vc_list1
            good_list = [good_list,n];
        end
    catch
    end
    
    
end

%% 


% init variables
X = {};
X2 = {};
vc_ind = [];





for p = 1:length(vc_list1)
    X{p} = [];
    X2{p} = [];
    vc_cf = vc_list1(p);
    
    for n = 1:length(good_list)
        % case for trills
%         vc_ind1 = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf);
%         vc_ind2 = find(SUrate{good_list(n)}{1}.stim(:,2)>vc_list2(1) & SUrate{good_list(n)}{1}.stim(:,2)<vc_list2(2));
%         vc_ind3 = find(SUrate{good_list(n)}{1}.stim(:,3)>vc_list3(1) & SUrate{good_list(n)}{1}.stim(:,3)<vc_list3(2));
%         
%         vc_ind = intersect(vc_ind1,vc_ind2);
%         vc_ind = intersect(vc_ind,vc_ind3);
        % case for phees
        vc_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == vc_cf);
        pre_stim = SUrate{good_list(n)}{1}.xb.pre_stimulus_record_time;
        post_stim = 500;% SUrate{good_list(n)}{1}.xb.post_stimulus_record_time;
        stim_dur = SUrate{good_list(n)}{1}.xb.stimulus_ch1(1,5);
        trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
        PSTH = [];
        for v = 1:length(vc_ind)
            PSTH = [PSTH; mean(SUrate{good_list(n)}{1}.PSTH{vc_ind(v)}(:,1:trial_dur),1)];
        end
        PSTH = mean(PSTH,1);
%         PSTH = PSTH/(range(PSTH)+5);
        X{p} = [X{p};PSTH];
        X2{p} = [X2{p};PSTH-mean(SUrate{good_list(n)}{1}.spont)];
    end
end


figure
cmap = colormap(jet(length(vc_list1)));
for p = 1:length(vc_list1)
        xs = 1:trial_dur;
        h = 30;
        for i = 1:trial_dur
            ys(i) = gaussian_kern_reg(xs(i),xs,mean(X2{p},1),h);
        end
        plot(xs-pre_stim,ys,'LineWidth',2,'Color',cmap(p,:))
        hold on
end
legend({'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0'});
ylabel('FR (spks/s)');
xlabel('time (ms)')


%% finding normalizing factor
norm_fact = zeros(1,length(good_list));
tic
for n = 1:length(good_list)
    figure(n)
    if mod(n,10) ==1
        fprintf(['%4d /' num2str(length(good_list)) ' time : %6.2f sec \n'],n,toc')
    end
    PSTH2 = [];
    for st = 1:length(vc_list1)
        xs = 1:length( X{st}(n,:));
        h = 20;
        ys = [];
        for i = 1:length( X{st}(n,:))
            ys(i) = gaussian_kern_reg(xs(i),xs,X{st}(n,:),h);
        end
        PSTH2 = [PSTH2; ys];
        if isempty(intersect(n,test_ind))
            plot(ys)
            hold on
        end
    end
    hold off
    drawnow
    norm_fact(n) = range(reshape(PSTH2,[],1));
    if norm_fact(n) < 1
        norm_fact(n) =1;
    end

end


for p = 1:length(vc_list1)
    for n = 1:length(good_list)
        X2{p}(n,:) = X2{p}(n,:)/(norm_fact(n)+5); 
    end
end

%% OLE and categorization
% 03/08/2021

%{ 


run PSTH first 
Training set SD 
[0 0.25 0.5] is 1
[4.5 4.75 5] is 0

%}

% N = length(SUrate);
% phee_list = round([7.16:0.12:9.56]*100)/100;
% pre_stim = 300;
% post_stim = 500;
% stim_dur = 1180;
% trial_dur = pre_stim+ post_stim + stim_dur; 


train_ind = randsample(10,5);
test_ind = setdiff(1:10,train_ind);
% set2 = [10.33,10.74];
% set1 = [6.64,7.05];
bothset = [set1 set2];
allset = 0:0.25:5;
C_train = [ones(length(set1),1); zeros(length(set2),1)].';

twindow = 100;
totalP = 100;

% W = zeros(trial_dur-twindow); %,N,totalP);
% C_test= zeros(trial_dur-twindow,length(bothset));% ,totalP);


lambda = 0.1; %L2 regularization constant
tic
W_all = {};
C_test_all = {};

%%

goodN = length(good_list);
parfor p = 1:totalP
    R = zeros(goodN,length(set1) + length(set2));
    R_test = zeros(goodN,length(set1) + length(set2));
    W = zeros(trial_dur-twindow,goodN);%,totalP);
    C_test= zeros(trial_dur-twindow,length(bothset));% ,totalP);
    %     fprintf(['%4d /' num2str(totalP) ' permutations.'],p) %  time : %6.2f sec \n'],p,toc')
    disp(p)
    train_ind = randsample(10,5);
    test_ind = setdiff(1:10,train_ind);
    for t = 1:trial_dur-twindow
        for n = 1:goodN
            for s = 1:length(bothset)
                vc_cf = bothset(s);
                % for trills
                %                 vc_ind1 = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf);
                %                 vc_ind2 = find(SUrate{good_list(n)}{1}.stim(:,2)>vc_list2(1) & SUrate{good_list(n)}{1}.stim(:,2)<vc_list2(2));
                %                 vc_ind3 = find(SUrate{good_list(n)}{1}.stim(:,3)>vc_list3(1) & SUrate{good_list(n)}{1}.stim(:,3)<vc_list3(2));
                %
                %                 vc_ind = intersect(vc_ind1,vc_ind2);
                %                 vc_ind = intersect(vc_ind,vc_ind3);
                %                 vc_ind = vc_ind(1); % change to incorporate all
                % for phees
                
                vc_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == vc_cf);
                R(n,s) =  mean2(SUrate{good_list(n)}{1}.PSTH{vc_ind}(train_ind,t:t+twindow)); %/(norm_fact(n)+5);
                R(n,s) =  mean2(SUrate{good_list(n)}{1}.PSTH{vc_ind}(train_ind,t:t+twindow)); %/(norm_fact(n)+5);
                R_test(n,s)= mean2(SUrate{good_list(n)}{1}.PSTH{vc_ind}(test_ind,t:t+twindow)); %/(norm_fact(n)+5));
            end
        end
        W_hat = C_train*((R.'*R + lambda^2)\R.');
        W(t,:) =  W_hat;
        C_test(t,:) =  W_hat*R_test;
        
    end
    W_all{p} = W;
    C_test_all{p} = C_test;
end
toc


W = zeros(trial_dur-twindow,goodN,totalP);
C_test= zeros(trial_dur-twindow,length(bothset),totalP);
for p = 1:totalP
    C_test(:,:,p) = C_test_all{p};
    W(:,:,p) = W_all{p};
end
C_test_mean = mean(C_test,3);

figure(1)
trial_time_x =[1:trial_dur-twindow]-pre_stim+twindow/2;
for s = 1:6
    plot(trial_time_x,C_test_mean(:,s))
    hold on
end
axis([-inf, inf, -inf, inf])
hold off

drawnow

%%
% test_set =  round([7.46:0.41:9.92]*100)/100;
% testing the other range with W

W_mean = mean(W,3);
C_test2_all = {};
parfor p = 1:totalP
    %         train_ind = randsample(10,5);
    %     test_ind = setdiff(1:10,train_ind);
    test_ind = randsample(10,7);
    disp(p)
    C_test2 = zeros(trial_dur-twindow,length(test_set));
    R_test2 = zeros(goodN,length(test_set));
    for t = 1:trial_dur-twindow
        for n = 1:goodN
            for s = 1:length(test_set)
                vc_cf = test_set(s);
                %                 vc_ind1 = find(SUrate{good_list(n)}{1}.stim(:,1) == vc_cf);
                %                 vc_ind2 = find(SUrate{good_list(n)}{1}.stim(:,2)>vc_list2(1) & SUrate{good_list(n)}{1}.stim(:,2)<vc_list2(2));
                %                 vc_ind3 = find(SUrate{good_list(n)}{1}.stim(:,3)>vc_list3(1) & SUrate{good_list(n)}{1}.stim(:,3)<vc_list3(2));
                %
                %                 vc_ind = intersect(vc_ind1,vc_ind2);
                %                 vc_ind = intersect(vc_ind,vc_ind3);
                
                vc_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == vc_cf);
                
                vc_ind = vc_ind(1); % change to incorporate all
                %             R(n,s) =  mean2(SUrate{n}{1}.PSTH{phee_ind}(train_ind,t:t+twindow));
                R_test2(n,s)= mean2(SUrate{good_list(n)}{1}.PSTH{vc_ind}(test_ind,t:t+twindow)); % /norm_fact(n);
            end
        end
        %     W_hat = C_train*((R.'*R + lambda^2)\R.');
        %     W(t,:,p) =  W_hat;
        C_test2(t,:) =  W_mean(t,:)*R_test2;
        
        
    end
    C_test2_all{p} = C_test2;
end

C_test2= zeros(trial_dur-twindow,length(test_set),totalP);
for p = 1:totalP
    C_test2(:,:,p) = C_test2_all{p};
end

C_test2_mean = mean(C_test2,3);

%% figures
cmap = colormap(jet(length(test_set)));
figure()
for s = 1:length(test_set)
    xs = 1:length(C_test2_mean(:,1));
    h = 10;
    ys = [];
    for i = 1:length(C_test2_mean(:,1))
        ys(i) = gaussian_kern_reg(xs(i),xs,C_test2_mean(:,s).',h);
    end
    plot(ys,'Color',cmap(s,:))
    hold on
    drawnow
end





figure
D_stim = mean(C_test2(pre_stim+600:pre_stim+stim_dur,:,:),1);
D_stim_mean = mean(D_stim,3);
D_stim_err = std(D_stim,[],3);%/sqrt(totalP);
errorbar(test_set,D_stim_mean,D_stim_err);
hold on
D_post = mean(C_test2(pre_stim+stim_dur:pre_stim+stim_dur+200,:,:),1);
D_post_mean = mean(D_post,3);
D_post_err = std(D_post,[],3);%/sqrt(totalP);

errorbar(test_set,D_post_mean,D_post_err);




% figure
% imagesc(C_test2_mean);
% sfig = surf(C_test2_mean);
% sfig.EdgeColor = 'none';
%% 
filt_im = imgaussfilt(C_test2_mean,[4,0.8]);
figure
cmap = colormap(jet(length(test_set)));


time_axis = [1:trial_dur-twindow]-pre_stim+twindow/2;
sfig = surf(test_set,time_axis,filt_im);
xlabel('Center Frequency, (kHz)');
ylabel('time');
sfig.EdgeColor = 'none';
imagesc(test_set,time_axis,filt_im);
caxis([0,1])


%% further testing

% W_mean_abs = abs(W_mean);
imagesc(W_mean);
test = mean(abs(W_mean),1);
bar(test)
histogram(test,14);
test_ind = find(test < 5*1e-4);
length(test_ind);



% end Optimal Linear Estimator









