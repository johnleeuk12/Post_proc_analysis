


%% PSTH
N = length(SUrate);
phee_list = round([7.16:0.12:9.56]*100)/100;

% good_list = [];
% for n = 1:N
%     counter = 0;
%     for p = 1: length(phee_list)
%         phee_cf = phee_list(p);
%         spont =  mean(SUrate{n}{1}.spont);       
%         SD = std(SUrate{n}{1}.spont,0,1)/sqrt(41);
%         
%         if SUrate{n}{1}.mean(p) > spont + 2*SD 
%             counter = counter+1;
%         end
%     end
%     if counter >1
%         good_list = [good_list;n];
%     end
% end


good_list = 1:N;
for p = 1:length(phee_list)
    X{p} = [];
    X2{p} = [];
    phee_cf = phee_list(p);
    
    for n = 1:length(good_list)
        phee_ind = find(SUrate{good_list(n)}{1}.xb.stimulus_ch1(:,10) == phee_cf);
        
        pre_stim = SUrate{good_list(n)}{1}.xb.pre_stimulus_record_time;
        post_stim = 500;
        stim_dur = 1180;
        %Gaussian smoothing
        trial_dur = pre_stim + stim_dur + post_stim;% + post_stim;
        PSTH = mean(SUrate{good_list(n)}{1}.PSTH{phee_ind}(:,1:trial_dur),1);
%         PSTH = PSTH/(range(PSTH)+5);
        X{p} = [X{p};PSTH];
        X2{p} = [X2{p};PSTH-mean(SUrate{good_list(n)}{1}.spont)];
    end
end

figure
cmap = colormap(jet(21));
for p = 1:2:length(phee_list)
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
norm_fact = zeros(1,N);
tic
for n = 1:N
    if mod(n,10) ==1
        fprintf(['%4d /' num2str(N) ' time : %6.2f sec \n'],n,toc')
    end
    PSTH2 = [];
    for st = 1:21
        xs = 1:length( X{st}(n,:));
        h = 10;
        ys = [];
        for i = 1:length( X{st}(n,:))
            ys(i) = gaussian_kern_reg(xs(i),xs,X{st}(n,:),h);
        end
        PSTH2 = [PSTH2; ys];
        norm_fact(n) = range(reshape(PSTH2,[],1));
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

N = length(SUrate);
phee_list = round([7.16:0.12:9.56]*100)/100;
pre_stim = 300;
post_stim = 500;
stim_dur = 1180;
trial_dur = pre_stim+ post_stim + stim_dur; 


train_ind = randsample(10,5);
test_ind = setdiff(1:10,train_ind);
set1 = [7.16,7.28,7.40];
set2 = [9.32,9.44,9.56];
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
parfor p = 1:totalP
    R = zeros(N,length(set1) + length(set2));
    R_test = zeros(N,length(set1) + length(set2));
    W = zeros(trial_dur-twindow,N);%,totalP);
    C_test= zeros(trial_dur-twindow,length(bothset));% ,totalP);
    %     fprintf(['%4d /' num2str(totalP) ' permutations.'],p) %  time : %6.2f sec \n'],p,toc')
    disp(p)
    train_ind = randsample(10,5);
    test_ind = setdiff(1:10,train_ind);
    for t = 1:trial_dur-twindow
        for n = 1:N
            for s = 1:length(bothset)
                phee_cf = bothset(s);
                phee_ind = find(SUrate{n}{1}.xb.stimulus_ch1(:,10) == phee_cf);
                R(n,s) =  mean2(SUrate{n}{1}.PSTH{phee_ind}(train_ind,t:t+twindow))/(norm_fact(n)+5);
                R_test(n,s)= mean2(SUrate{n}{1}.PSTH{phee_ind}(test_ind,t:t+twindow)/(norm_fact(n)+5));
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

%%
W = zeros(trial_dur-twindow,N,totalP);
C_test= zeros(trial_dur-twindow,length(bothset),totalP);
for p = 1:totalP
    C_test(:,:,p) = C_test_all{p};
    W(:,:,p) = W_all{p};
end
C_test_mean = mean(C_test,3);

figure(1)
for s = 1:6
    plot(C_test_mean(:,s))
    hold on
end

drawnow
%%
test_set =  round([7.52:0.12:9.2]*100)/100;
% testing the other range with W

W_mean = mean(W,3);
C_test2_all = {};
parfor p = 1:totalP
    %         train_ind = randsample(10,5);
    %     test_ind = setdiff(1:10,train_ind);
    test_ind = randsample(10,7);
    disp(p)
    C_test2 = zeros(trial_dur-twindow,length(test_set));
    R_test2 = zeros(N,length(test_set));
    for t = 1:trial_dur-twindow
        for n = 1:N
            for s = 1:length(test_set)
                phee_cf = test_set(s);
                phee_ind = find(SUrate{n}{1}.xb.stimulus_ch1(:,10) == phee_cf);
                %             R(n,s) =  mean2(SUrate{n}{1}.PSTH{phee_ind}(train_ind,t:t+twindow));
                R_test2(n,s)= mean2(SUrate{n}{1}.PSTH{phee_ind}(test_ind,t:t+twindow))/norm_fact(n);
            end
        end
        %     W_hat = C_train*((R.'*R + lambda^2)\R.');
        %     W(t,:,p) =  W_hat;
        C_test2(t,:) =  W_mean(t,:)*R_test2;
        
        
    end
    C_test2_all{p} = C_test2;
end
%%
C_test2= zeros(trial_dur-twindow,length(test_set),totalP);

for p = 1:totalP
        C_test2(:,:,p) = C_test2_all{p};
end

figure
C_test2_mean = mean(C_test2,3);
D_stim = mean(C_test2(pre_stim+0:pre_stim+300,:,:),1);
D_stim_mean = mean(D_stim,3);
D_stim_err = std(D_stim,[],3);%/sqrt(totalP);
errorbar(test_set,D_stim_mean,D_stim_err);
hold on
D_post = mean(C_test2(pre_stim+stim_dur:trial_dur-twindow,:,:),1);
D_post_mean = mean(D_post,3);
D_post_err = std(D_post,[],3);%/sqrt(totalP);

errorbar(test_set,D_post_mean,D_post_err);


cmap = colormap(jet(15));
figure()
for s = 1:15
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

% figure
% imagesc(C_test2_mean);
% sfig = surf(C_test2_mean);
% sfig.EdgeColor = 'none';

filt_im = imgaussfilt(C_test2_mean,[4,0.8]);
figure


time_axis = [1:trial_dur-twindow]-pre_stim+twindow/2;
sfig = surf(test_set,time_axis,filt_im);
xlabel('Center Frequency, (kHz)');
ylabel('time');
sfig.EdgeColor = 'none';
imagesc(test_set,time_axis,filt_im);
caxis([0,1])



% end Optimal Linear Estimator


%% Neural dynamics
N = length(SUrate);
figure;
set(gcf, 'Position', [400 200 800 800]);

phee_list = round([7.16:0.24:9.56]*100)/100;
X = {};
tic
for p = 1:length(phee_list)
    phee_cf = phee_list(p);
    X{p} = [];
    fprintf(['%4d /' num2str(length(phee_list)) ' time : %6.2f sec \n'],p,toc)
    
    for n = 1:N
        %         if mod(n,10) ==1
        %             fprintf(['%4d /' num2str(N) ' time : %6.2f sec \n'],n,toc')
        %         end
        
        phee_ind = find(SUrate{n}{1}.xb.stimulus_ch1(:,10) == phee_cf);
        pre_stim = SUrate{n}{1}.xb.pre_stimulus_record_time;
        post_stim = 500;
        stim_dur = 1180;
        %Gaussian smoothing
        trial_dur = pre_stim + stim_dur;% + post_stim;
        xs = 1:trial_dur;
        h = 30;
        for i = 1:trial_dur
            ys(i) = gaussian_kern_reg(xs(i),xs,mean(SUrate{n}{1}.PSTH{phee_ind}(:,1:trial_dur),1),h);
        end
        %             ys = movmean(mean(SUrate{n}{1}.PSTH{phee_ind}(:,1:trial_dur),1),100);
        ys = ys-mean(SUrate{n}{1}.spont);
        ys = ys/(range(ys)+5);
        X{p} = [X{p}; ys];
    end
    

end
% colormap default
% patch(T1,T2,T3);
%%
figure
pl = {};
cmap = colormap(jet(11));
for p = 1:2:11
        coeff = pca(X{p});
    T1 = coeff(:,1);
    T2 = coeff(:,2);
    T3 = coeff(:,3);
    
    pl{p} = plot3(T1,T2,T3,'Color',cmap(p,:));
    hold on
    scatter3(T1(1),T2(1),T3(1),100,cmap(p,:),'*')
    scatter3(T1(end),T2(end),T3(end),100,cmap(p,:),'o','filled')
%     scatter3(T1(pre_stim),T2(pre_stim),T3(pre_stim),200,cmap(p,:),'^','filled')
%     scatter3(T1(pre_stim+stim_dur),T2(pre_stim+stim_dur),T3(pre_stim+stim_dur),200,cmap(p,:),'v','filled')
    drawnow
    pause(0.1)
    
end

legend([pl{1} pl{3} pl{5} pl{7} pl{9} pl{11}])