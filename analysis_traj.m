%% PSTH from DD 

nrep = 10;
PSTH = zeros(21,1980);
for st = 1:21
    cumu_data = [];
    for r = 1:10
    cumu_data = [cumu_data; mean(DD(r+(st-1)*10).data,1)];
    end
    PSTH(st,:) = smoothdata(mean(cumu_data,1),'gaussian',100);
end
cmap = parula(21);
figure
for st = 1:2:21
    plot(PSTH(st,:),'Color',cmap(st,:))
    hold on
end
hold off

%% calculating distance of different PCs from 0SD
X = {};
figure
cmap = parula(21);
for pc = 1:8
    subplot(2,5,pc)
    for st = 1:21
        plot(1:99,Dtraj(st).data(pc,:),'Color',cmap(st,:));
        hold on
    end
    xline(15);
    xline(74);
    axis([0, 100,min(min(Dtraj(1).data(pc,:))*2,-1.5),max(max(Dtraj(1).data(pc,:))*2,1.5)])
    title(['pc' num2str(pc)])
    hold off
end

%%

% stim_set = 7.28:0.12:9.56;
stim_set = 4.76:0.12:7.04;
for pc = 1:8
    X{pc}.data = zeros(20,99);
    for st = 1:20
        for t = 1:99
            X{pc}.data(st,t) = abs(Dtraj(21).data(pc,t)-Dtraj(st).data(pc,t));
        end
    end
    X{pc}.data = X{pc}.data/(max(max(X{pc}.data))-min(min(X{pc}.data)));
    figure(pc)
    subplot(1,2,1)
    imagesc(stim_set,linspace(-300,1680,99),X{pc}.data.')
    title(num2str(pc))
    colorbar
    drawnow
end

%%

% plotting distance at different times for each PC

% x_axis = 7.28:0.12:9.56;
x_axis = 4.76:0.12:7.04;
for pc = 1:8
    %     subplot(2,4,pc)
    figure(pc)
    subplot(1,2,2)
    X{pc}.onset = mean(X{pc}.data(:,15:24),2);
    X{pc}.mid = mean(X{pc}.data(:,35:54),2);
    X{pc}.end = mean(X{pc}.data(:,55:74),2);
    X{pc}.offset = mean(X{pc}.data(:,80:90),2);
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


