%% Population Analysis .

N_pool = size(Pool,2);
for ss = 1:N_pool
    N_list{ss} = unique([Pool{ss}.neuron_nb]);
end

% find intersecting list 
N_list_int = N_list{1};
if size(Pool,2) >1
    for ss = 1:N_pool-1
        N_list_int = intersect(N_list_int,N_list{ss+1});
    end
end
N_list_int = N_list_int(2:end);
% SUrate{1, 1}{1, 1}.mean
% nat_ind = [1,3,5,8,10,11,13,15,17,19]; % for twitters only. seems like stim 4 and 5 were reversed
nat_ind = 1:2:19;

ind_1 = [];
ind_2 = [];
cat_1 = {};
cat_1.pool = [];
cat_2 = {};
cat_2.pool = [];

remain = [];
ss =3;
for n = 1:length(N_list_int)
%     if length(SUrate{n}{ss}.mean) ~=16 % this is only for twitter CF 
        vq =[];
        %     if length(SUrate{n}{ss}.mean) == 21 || length(SUrate{n}{ss}.mean) == 16
        
        x = SUrate{n}{ss}.mean-mean(SUrate{n}{ss}.spont);
        if length(SUrate{n}{ss}.mean) == 21
            %     elseif length(SUrate{n}{ss}.mean) == 16
            %         x = [x;nan;nan;nan;nan;nan];
        elseif length(SUrate{n}{ss}.mean) == 23
            %         vq2 = interp1(x,v,xq,'spline');
            vq = interp1(1:2:21,x(1:11),1:21,'spline');
            x = vq.';
        end
        
        x = x/max(abs(x));
        if SUrate{n}{2}.mean(nat_ind)> mean(SUrate{n}{2}.spont)% + std(SUrate{n}{2}.spont)/sqrt(20)
            ind_1 = [ind_1 n];
            cat_1.pool = [cat_1.pool x];
        else
            ind_2 = [ind_2 n];
            cat_2.pool = [cat_2.pool x];
        end
%     end
end

%% Clustering for tuning

%calculate euc distance between each tuning curves
distmap = zeros(size(cat_1.pool,2),size(cat_1.pool,2));
for f1 = 1:size(cat_1.pool,2)
    for f2 = 1:size(cat_1.pool,2)
        if f1 ~= f2
            distmap(f1,f2) = pdist([cat_1.pool(:,f1),cat_1.pool(:,f2)].');
        else
            distmap(f1,f2) = 0;
        end
    end
end

% re-ordering using hierarchical clustering
Z = linkage(cat_1.pool.','average');
leafOrder = optimalleaforder(Z,distmap);

new_pool = [];
for f = 1:length(leafOrder)
    new_pool = [new_pool, cat_1.pool(:,leafOrder(f))];
end


new_pool = new_pool(:,1:end-1);
if ss ==2
    stim_index = Pool{2}(1).xb.stimulus_ch1(:,12);
else
    stim_index = 4:0.4:12;
    
end
figure
imagesc(stim_index,[],new_pool.')
%%

% Cat 2

distmap = zeros(size((cat_2.pool),2),size((cat_2.pool),2));
for f1 = 1:size((cat_2.pool),2)
    for f2 = 1:size((cat_2.pool),2)
        if f1 ~= f2
            distmap(f1,f2) = pdist([cat_2.pool(:,f1),cat_2.pool(:,f2)].');
        else
            distmap(f1,f2) = 0;
        end
    end
end

% re-ordering using hierarchical clustering
Z = linkage(cat_2.pool.','average');
leafOrder2 = optimalleaforder(Z,distmap);

new_pool2 = [];
for f = 1:length(leafOrder2)
    new_pool2 = [new_pool2, cat_2.pool(:,leafOrder2(f))];
end


new_pool2 = new_pool2(:,1:end-1);
if ss ==2
    stim_index = Pool{2}(1).xb.stimulus_ch1(:,12);
else
    stim_index = 4:0.4:12;
    
end
figure
imagesc(stim_index,[],new_pool2.')




%%
% population average
belt.raw = new_pool;
belt.mean = mean(new_pool(:,1:42),2);
belt.err = std(new_pool(:,1:42),[],2);

% core.raw = new_pool;
% core.mean = mean(new_pool(:,1:42),2);
% core.err = std(new_pool(:,1:42),[],2);

% plot(stim_index,cat_1.mean)
figure
errorbar(stim_index,core.mean,core.err/sqrt(length(core.err)))
hold on
errorbar(stim_index,belt.mean,belt.err/sqrt(length(belt.err)))
% figure
% bar(stim_index,cat_1.err)






%% find units that respond to real vocalizations

% nat_ind = [1,3,5,8,10,11,13,15,17,19]; % for twitters only. seems like stim 4 and 5 were reversed
nat_ind = 1:2:19;

ind_1 = [];
ind_2 = [];


remain = [];
ss =3;
for n = 1:length(SUrate)
    count_sig = 0;
    for id = 1:length(nat_ind)
        if SUrate{n}{2}.mean(id)> mean(SUrate{n}{2}.spont)% + std(SUrate{n}{2}.spont)/sqrt(20)
            count_sig = count_sig +1;
        end
    end
    if count_sig >1
        ind_1 = [ind_1 n];
    else
        ind_2 = [ind_2 n];
    end
    
end


%
%% Analysis onset, sustained and offset responses 
close all

PreStim = 300;
PostStim = 700;
Stimdur = 2180-PostStim-PreStim;
N_stim = 16; %16 for phee10 21 for Twitter and phee12
ss = 3;
ana = {};
ana.onset =  [];
ana.sust = [];
ana.offset = [];
for nu = 1:length(ind_1)
    
    
    for st = 1:N_stim
        ana.onset.mean(nu,st) = mean2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim:PreStim+250));
        ana.onset.error(nu,st) = std2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim:PreStim+250))/sqrt(250*16);
        
        ana.sust.mean(nu,st) = mean2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim+250:PreStim+Stimdur));
        ana.sust.error(nu,st) = std2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim+250:PreStim+Stimdur))/sqrt(850*16);
        
        ana.offset.mean(nu,st) = mean2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim+Stimdur:PreStim+Stimdur+250));
        ana.offset.error(nu,st) = std2(SUrate{ind_1(nu)}{ss}.PSTH{st}(:,PreStim+Stimdur:PreStim+Stimdur+250))/sqrt(250*16);
    end
end


PreStim = 200;
Stimdur = 100;
PostStim = 300;

% tuning to bp noise and to PT 
for nu = 1:length(ind_1)
    for st = 1:length(SUrate{1}{1}.mean)
        ana.PT.mean(nu,st) = mean2(SUrate{ind_1(nu)}{1}.PSTH{st}(:,PreStim:PreStim+Stimdur+100)); %getting also offset responses
        ana.PT.error(nu,st) = std2(SUrate{ind_1(nu)}{1}.PSTH{st}(:,PreStim:PreStim+Stimdur+100))/sqrt(size(SUrate,2)*(PreStim+Stimdur+100));
        ana.BP.mean(nu,st) = mean2(SUrate{ind_1(nu)}{4}.PSTH{st}(:,PreStim:PreStim+Stimdur+100)); %getting also offset responses
        ana.BP.error(nu,st) = std2(SUrate{ind_1(nu)}{4}.PSTH{st}(:,PreStim:PreStim+Stimdur+100))/sqrt(size(SUrate,2)*(PreStim+Stimdur+100));
    end
end

stim_label1 = 4:0.4:10; %10 for phee for hole 15, 12 for others
stim_labelPT = logspace(log10(4),log10(32),31);
for nu = 1:length(ind_1)
    figure(nu)
    errorbar(stim_label1,ana.onset.mean(nu,:),ana.onset.error(nu,:),'LineWidth',2);
    hold on
    errorbar(stim_label1,ana.sust.mean(nu,:),ana.sust.error(nu,:),'LineWidth',2);
    errorbar(stim_label1,ana.offset.mean(nu,:),ana.offset.error(nu,:),'LineWidth',2);
    % BP and PT
    errorbar(stim_labelPT,ana.PT.mean(nu,:),ana.PT.error(nu,:),':','LineWidth',2);
    errorbar(stim_labelPT,ana.BP.mean(nu,:),ana.BP.error(nu,:),':','LineWidth',2);
    legend({'onset','sustained','offset','PT','BP'});
    axis([4 10 -inf inf])
    title(['unit' num2str(SUrate{ind_1(nu)}{1}.nid)])
    xlabel('CF(kHz)');
    ylabel('FR (spk/s)');
    drawnow()
end



%% hierarchical clustering

for nu = 1:length(ind_1)
    M.onset(nu,:) = ana.onset.mean(nu,:)/max(ana.onset.mean(nu,:));
    M.sust(nu,:) = ana.sust.mean(nu,:)/max(ana.sust.mean(nu,:));
    M.offset(nu,:) = ana.offset.mean(nu,:)/max(ana.offset.mean(nu,:));
    M.PT(nu,:) = ana.PT.mean(nu,:)/max(ana.PT.mean(nu,1:15)); %11:25
    M.BP(nu,:) = ana.BP.mean(nu,:)/max(ana.BP.mean(nu,1:15));
end

% normalized population average
% figure
% plot(stim_label1,mean(M.offset,1,'omitnan'),'LineWidth',2)
% hold on
% plot(stim_label1,mean(M.sust,1,'omitnan'),'LineWidth',2)
% plot(stim_label1,mean(M.onset,1,'omitnan'),'LineWidth',2)
% plot(stim_labelPT,mean(M.PT,1,'omitnan'),':','LineWidth',2)
% plot(stim_labelPT,mean(M.BP,1,'omitnan'),':','LineWidth',2)
%     legend({'onset','sustained','offset','PT','BP'});
%     axis([4 10 -inf inf])

% plot(stim_label1,std(M.offset,1,'omitnan')/sqrt(size(SUrate,2)));
T = M.onset;%([1:4,6:16],:);
% T = M.offset([1:3,5:end],:);
T2 = M.PT(:,1:15);
% T = M.BP(:,11:25); %(:,11:25);
N = size(T,1);


distmap = zeros(N,N);
for f1 = 1:N
    for f2 = 1:N
        if f1 ~= f2
            distmap(f1,f2) = pdist([T(f1,:);T(f2,:)]);
        else
            distmap(f1,f2) = 0;
        end
    end
end

% re-ordering using hierarchical clustering
Z = linkage(T,'average');
leafOrder2 = optimalleaforder(Z,distmap);

new_pool2 = [];
new_pool3 = [];
for f = 1:length(leafOrder2)
    new_pool2 = [new_pool2; T(leafOrder2(f),:)];
    new_pool3 = [new_pool3; T2(leafOrder2(f),:)];
end


% new_pool2 = new_pool2(:,1:end-1);
% if ss ==2
%     stim_index = Pool{2}(1).xb.stimulus_ch1(:,12);
% else
%     stim_index = 4:0.4:12;
%     
% end
% figure
% imagesc(stim_label,[],T);
figure
imagesc(stim_label1,[],new_pool2)
xlabel('CF (kHz)')
ylabel('units')

figure
imagesc(stim_label1,[],new_pool3)
xlabel('CF (kHz)')
ylabel('units')


%% Analysis VT Phee
N = length(SUrate);
phee_cf = 7.16;
X = [];
for n = 1:N
    max_ind = find(SUrate{n}{1}.mean == max(SUrate{n}{1}.mean));
    max_ind = max_ind(1);
    spont = mean(SUrate{n}{1}.spont);
    SD = std(SUrate{n}{1}.spont,0,1);
    if SUrate{n}{1}.mean(max_ind)>2
        if SUrate{n}{1}.mean(max_ind) > spont +2*SD
            if max(SUrate{n}{2}.mean) > 2
                phee_ind = find(SUrate{n}{1, 2}.xb.stimulus_ch1(:,10) == phee_cf);
                if ~isempty(phee_ind)
                    X = [X; [SUrate{n}{1}.xb.stimulus_ch1(max_ind,8)*1e-3 (SUrate{n}{2}.mean(phee_ind)-mean(SUrate{n}{2}.spont))/max(SUrate{n}{2}.mean) SUrate{n}{2}.nid]];
                end
            end
        end
    end
end

[B,I] = sort(X(:,1));
X2 = [B,X(I,2)];
M = movmedian(X2(:,2),10);
figure
hold off
scatter(X2(:,1),X2(:,2))
hold on
plot(B,M);
xline(phee_cf)
set(gca, 'xScale', 'log')
xticks([0 : 3.5:30])
axis([2 inf -1 1])
title([num2str(phee_cf) 'kHz'])


%% Analysis VT trill

N = length(SUrate);
X = [];
% trill_cf = [9.8, 10];
trill_cf = [9.8, 10];
trill_fm = [27,30];
trill_fmd = [925 1000];
for n = 1:N
    if ~isempty(SUrate{n}{1})
        max_ind = find(SUrate{n}{1}.mean == max(SUrate{n}{1}.mean));
        max_ind = max_ind(1);
        spont = mean(SUrate{n}{1}.spont);
        SD = std(SUrate{n}{1}.spont,0,1);
        
        trill_ind_cf = find(SUrate{n}{2}.stim(:,1) >trill_cf(1) & SUrate{n}{2}.stim(:,1) <trill_cf(2)...
            & SUrate{n}{2}.stim(:,3)> trill_fmd(1) & SUrate{n}{2}.stim(:,3) < trill_fmd(2));
        trill_ind_fm = find(SUrate{n}{2}.stim(:,2) > trill_fm(1) & SUrate{n}{2}.stim(:,2) < trill_fm(2)  ...
            & SUrate{n}{2}.stim(:,3)> trill_fmd(1) & SUrate{n}{2}.stim(:,3) < trill_fmd(2));
        trill_ind = intersect(trill_ind_cf, trill_ind_fm);
        if length(trill_ind) >1
            trill_ind = trill_ind(end);
        end
        if ~ isempty(trill_ind)
            if SUrate{n}{1}.mean(max_ind) > spont +2*SD
                if max(SUrate{n}{2}.mean(trill_ind_fm))>1
                    X = [X; [SUrate{n}{1}.xb.stimulus_ch1(max_ind,8)*1e-3 (SUrate{n}{2}.mean(trill_ind)-mean(SUrate{n}{2}.spont(trill_ind_fm)))/max(SUrate{n}{2}.mean(trill_ind_fm))]];
                end
            end
        end
    end
    
end


figure
[B,I] = sort(X(:,1));
X2 = [B,X(I,2)];
M = movmean(X2(:,2),5);
hold off
scatter(X2(:,1),X2(:,2))
hold on
% xline(6.64);
xline(mean(trill_cf));
plot(B,M);

set(gca, 'xScale', 'log')
xticks([0 : 3.5:30])
title([num2str(mean(trill_cf)) 'kHz'])

