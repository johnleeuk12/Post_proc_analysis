%% reformatting data for trills

%% VT type

VT_code = 31;
switch VT_code
    case 11
        trial_dur = 1980;
        stim_set = 7.28:0.12:9.56;
        stim_num = 21;
        
    case 12
        trial_dur = 1980;
        stim_set = 4.76:0.12:7.04;
        stim_num = 21;
        
    case 21
        trial_dur = 1206;
        stim_num = 11;
        stim_set = 6.64:0.41:10.74;
        newDtraj = {};
        
        newDtraj = Dtraj(3:end);
        newDtraj(10:11) = Dtraj(1:2);
        Dtraj = {};
        Dtraj = newDtraj;
    case 22
    case 31
        Dtraj_old = Dtraj;
        Dtraj(2:end) = Dtraj(1:end-1);
        Dtraj(1) = Dtraj_old(11);
        stim_num = 11;
        stim_set = 9.72:0.735:17.07;
%         stim_set = 2.37:0.735:9.72;
        

end

%%


nrep = 8;
PSTH = zeros(stim_num,trial_dur);
for st = 1:stim_num
    cumu_data = [];
    for r = 1:nrep
    cumu_data = [cumu_data; mean(DD(r+(st-1)*nrep).data,1)];
    end
    PSTH(st,:) = smoothdata(mean(cumu_data,1),'gaussian',100);
end
cmap = parula(stim_num);
figure
for st = 1:stim_num
    plot(PSTH(st,:),'Color',cmap(st,:))
    hold on
end
hold off
%% 
X = {};
figure
cmap = parula(stim_num);
for pc = 1:10
    subplot(2,5,pc)
    for st = 1:stim_num
        plot(1:length(Dtraj(st).data),smoothdata(Dtraj(st).data(pc,:),'gaussian',1),'Color',cmap(st,:));
        hold on
    end
%     xline(15);
%     xline(74);
    axis([0, length(Dtraj(st).data),min(min(Dtraj(1).data(pc,:))*2,-1.5),max(max(Dtraj(1).data(pc,:))*2,1.5)])
    title(['pc' num2str(pc)])
    hold off
end

% 
%% calculating distance of different PCs from 0SD

for pc = 1:10
    X{pc}.data = zeros(stim_num-1,length(Dtraj(st).data));
    for st = 1:stim_num-1
        for t = 1:length(Dtraj(st).data)
            X{pc}.data(st,t) = abs(Dtraj(1).data(pc,t)-Dtraj(st+1).data(pc,t)); 
        end
    end
    X{pc}.data = X{pc}.data/(max(max(X{pc}.data))-min(min(X{pc}.data)));
    fi = figure(pc);
    set(fi, 'Position', [500 400 1500 800]);
    subplot(1,2,1)
    imagesc(stim_set,linspace(-300,trial_dur-300,length(Dtraj(st).data)),X{pc}.data.')
    title(num2str(pc))
    colorbar
    drawnow
end

distan = [];
% for pc = 1:10
%     distan = mean(X{pc}.data,2);

%%

% plotting distance at different times for each PC

% x_axis = 7.28:0.12:9.56;
% x_axis = 4.76:0.12:7.04;
x_axis = stim_set(2:end);

for pc = 5
    %     subplot(2,4,pc)
    figure(pc);
    subplot(1,2,2)
    X{pc}.onset = mean(X{pc}.data(:,15:24),2);
    X{pc}.mid = mean(X{pc}.data(:,35:54),2);
    X{pc}.end = mean(X{pc}.data(:,55:74),2);
    X{pc}.offset = mean(X{pc}.data(:,80:85),2);

%     X{pc}.onset = mean(X{pc}.data(:,16:22),2);
%     X{pc}.mid = mean(X{pc}.data(:,20:30),2);
%     X{pc}.end = mean(X{pc}.data(:,25:35),2);
%     X{pc}.offset = mean(X{pc}.data(:,38:55),2);
    plot(x_axis, X{pc}.onset)
    hold on
    plot(x_axis, X{pc}.mid)
    plot(x_axis, X{pc}.end)
    plot(x_axis, X{pc}.offset)
    title(num2str(pc))
    axis([-inf,inf,0,1])
    hold off
    drawnow
    legend({'onset','mid','end','offset'})
end


%%

