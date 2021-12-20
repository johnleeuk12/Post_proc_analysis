%% reformatting data for trills

%% VT type
PreStim = 300;
PostStim = 500;
nreps = 10;



VT_code = 10;
switch VT_code
    case 10
        stim_dur = 1180;
        stim_set = 4.76:0.12:9.56;
        stim_num = 41;
        zero_ind = 21;  
    
    case {01, 11} 
        stim_dur = 1180;
        stim_set = 6.68:0.12:9.56;
        stim_num = 25;
        zero_ind = 5;
        
    case {02, 12}
        stim_dur = 1180;
        stim_set = 4.76:0.12:7.64;
        stim_num = 25;
        zero_ind = 21;
        
    case 21
        %         trial_dur = 1206;
        stim_dur = 406;
        stim_num = 13;
        stim_set = 5.82:0.41:10.74;
        zero_ind = 3;
        %         newDtraj = {};
        
        %         newDtraj = Dtraj(3:end);
%         newDtraj(10:11) = Dtraj(1:2);
%         Dtraj = {};
%         Dtraj = newDtraj;
    case 22
        stim_dur = 406;
        stim_num = 13;
        stim_set = 2.54:0.41:7.46;
        zero_ind = 11;
    case 31
        %         Dtraj_old = Dtraj;
        %         Dtraj(2:end) = Dtraj(1:end-1);
        %         Dtraj(1) = Dtraj_old(11);
        stim_dur = 1161;
        trial_dur = stim_dur+PostStim+ PreStim;
        
        stim_num = 13;
        stim_set = 8.25:0.735:17.07;
        %         stim_set = 2.37:0.735:9.72;
        zero_ind = 3;
        
    case 32
        stim_num = 13;
                stim_dur = 1161;

        stim_set = 2.37:0.735:11.19;
        zero_ind = 11;
        
end
trial_dur = stim_dur+PostStim+ PreStim;

%%


nrep = nreps;
PSTH = zeros(stim_num,trial_dur);

keep_neurons_ind = find(keep_neurons == 1);

for st = 1:stim_num
    cumu_data = [];
    for r = 1:nrep
    cumu_data = [cumu_data; mean(DD(r+(st-1)*nrep).data(keep_neurons_ind,:),1)];
    end
    PSTH(st,:) = smoothdata(mean(cumu_data,1),'gaussian',200);
end
cmap = parula(stim_num);
% cmap = flip(parula(stim_num));
figure(53)
for st = 1:2:stim_num
    plot([1:trial_dur]-PreStim,PSTH(st,:)*1e3,'Color',cmap(st,:),'LineWidth',2)
    hold on
end
axis([-PreStim, trial_dur-PreStim, -inf, inf])
hold off
%% 
% X = {};
figure(30)
% cmap = flip(parula(stim_num));
cmap = parula(stim_num);
for pc = 1:5
    subplot(1,5,pc)
    for st = 1:stim_num
        plot(1:length(Dtraj(st).data),smoothdata(Dtraj(st).data(pc,:),'gaussian',15),'Color',cmap(st,:));
        hold on
    end
%     xline(15);
%     xline(74);
    axis([0, length(Dtraj(st).data),min(min(Dtraj(pc).data(pc,:))*2,-1.5),max(max(Dtraj(pc).data(pc,:))*2,1.5)])
    title(['pc' num2str(pc)])
    hold off
end

% 
%% calculating distance of different PCs from 0SD
% else_ind = 1:stim_num;
% else_ind(zero_ind) = [];
X = {};
for pc = 1:3
    X{pc}.data = zeros(stim_num,length(Dtraj(st).data));
    for st = 1:stim_num
        for t = 1:length(Dtraj(st).data)
            X{pc}.data(st,t) = abs(Dtraj(zero_ind).data(pc,t)-Dtraj(st).data(pc,t)); 
        end
    end
    X{pc}.data = X{pc}.data/(max(max(X{pc}.data))-min(min(X{pc}.data)));
    fi = figure(pc);
    set(fi, 'Position', [500 400 1500 600]);
    subplot(1,2,1)
    imagesc(stim_set,linspace(-300,trial_dur-300,length(Dtraj(st).data)),X{pc}.data.')
    colormap(parula)
    title(num2str(pc))
    colorbar
    drawnow
end

distan = [];
% for pc = 1:10
%     distan = mean(X{pc}.data,2);

%%
% determining epoch times
epoch_starts = ceil(Dtraj(1).epochStarts/20);


% plotting distance at different times for each PC


x_axis = stim_set;

for pc = 1:3
    %     subplot(2,4,pc)
    figure(pc);
    subplot(1,2,2)
    X{pc}.pre = mean(mean(X{pc}.data(:,1:epoch_starts(2)-1),2));
    X{pc}.onset = mean(X{pc}.data(:,epoch_starts(2)+5:epoch_starts(3)),2);
    %     X{pc}.mid = mean(X{pc}.data(:,epoch_starts(2)+10:epoch_starts(3)-5),2);
    %     X{pc}.end = mean(X{pc}.data(:,epoch_starts(3):epoch_starts(3)+10),2);
    X{pc}.offset = mean(X{pc}.data(:,epoch_starts(3)+5:epoch_starts(3)+15),2);
    
    
    plot(x_axis, X{pc}.onset,'-g','LineWidth',2)
    hold on
    %     plot(x_axis, X{pc}.mid)
    %     plot(x_axis, X{pc}.end)
    plot(x_axis, X{pc}.offset,'-r','LineWidth',2)
    yline(X{pc}.pre,'--k');
    
    title(num2str(pc))
    axis([-inf,inf,0,1])
    hold off
    drawnow
    %     legend({'onset','mid','end','offset'})
    legend('onset+sustained','offset')
end


%% evaluation on units, variance explained, plotting PSTH for individual units
figure(51)
lat_T = explained(1:18);
lat_T = cumsum(lat_T/sum(lat_T));
plot(lat_T,'*k')
hold on
yline(0.8)


%%


sel_pc = 3;

figure(54)
imagesc(C)
lat_C = C(:,sel_pc);
% plot(sort(lat_C));
lat_C = abs(lat_C);

[sorted_C, I] = sort(lat_C,'descend');
norm_C = cumsum(sorted_C)/sum(sorted_C);
plot(norm_C);
hold on
yline(0.8);
yline(0.5);


nid_list = good_list_nid(keep_neurons);

Sel_nidPC = nid_list(I(find(norm_C<0.8)));


%%


NN = length(SUrate);
tpoint1 = PreStim + stim_dur+ 100;
tpoint2 = tpoint1 + 400 ;
tdur = tpoint2-tpoint1;


% for nn = Sel_nidPC
for nn = 1:length(Sel_nid_all)
    %     nn = Sel_nid(21);
    
    figure(50)
    for n = 1:NN
        if SUrate{n}{1}.nid == Sel_nid_all(nn)
            
            DR = zeros(length(SUrate{n}{1}.mean),tdur+1);
            for st = 1:length(SUrate{n}{1}.mean)
                DR(st,:) = mean(SUrate{n}{1}.PSTH{st}(:,tpoint1:tpoint2));
            end
            
            
            nPSTH = zeros(length(SUrate{n}{1}.mean),trial_dur);
            for st = 1:length(SUrate{n}{1}.mean)
                nPSTH(st,:) = mean(SUrate{n}{1}.PSTH{st}(:,1:trial_dur));
            end
            
            PSTH2 = imgaussfilt(nPSTH,[3,20],'Padding',0);
%             subplot(2,2,1)
            imagesc([1:trial_dur]-PreStim,SUrate{n}{1}.stim(:,1),PSTH2);

            colormap(flipud(gray))
            xline(0);
            xline(stim_dur);
            ylabel('kHz');
            xlabel('ms');
            title(['M60Fu' num2str(Sel_nid_all(nn))]);
%             subplot(2,2,2)
%             DR_2 = mean(DR(17:end,:),2)-mean(SUrate{n}{1}.spont);
%             errorbar(stim_set,DR_2,std(DR(17:end,:),0,2)/sqrt(300))
%             if DR_2(5)>1*std(SUrate{n}{1}.spont)
%                 distance = abs(DR_2-DR_2(5));
%                 subplot(2,2,3)
%                 plot(stim_set,distance);
%             end

        end
    end
    pause
    clf
end



% for nn = 1:NN
%     if length(SUrate{nn}{1}.pid)<2
%         disp(nn)
%         pause
%     end
% end
% 
% 

%% 
% NN = length(SUrate);
% 
% for nn = Sel_nid
%     for n = 1:NN
%         if SUrate{n}{1}.nid == nn
%             
%             nPSTH = zeros(length(SUrate{n}{1}.mean),tdur);
%             for st = 1:length(SUrate{n}{1}.mean)
%                 nPSTH(st,:) = mean(SUrate{n}{1}.PSTH{st}(:,1:tdur));
%             end
%         end
%     end
% 
% 
% end
% 
% figure(50)
% plot(stim_set,mean(nPSTH(17:end,:),2))


















































