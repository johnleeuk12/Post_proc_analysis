%% Pool data

% clear all
Pool = {};
rate = {};
raster = {};


% [Pool{1}, rate{1}, raster{1}] = pool_data(2320,[1:5],40,60);
[Pool{1}, rate{1}, raster{1}] = pool_data(1,[2:3],60,1180);
% [Pool{2}, rate{2}, raster{2}] = pool_data(104,8,40);
% [Pool{1}, rate{1}, raster{1}] = pool_data(1,[10],40);
% [Pool{4}, rate{4}, raster{4}] = pool_data(62,8,40);
% [Pool{5}, rate{5}, raster{5}] = pool_data(2343,[12,13],40);

pause(0.1)

SUrate = {};
SUrate = get_tuning(Pool,rate,raster,0);
% use ana2Dstim for m60F, and ana2Dstimv2 for m160E
SUrate = ana_2Dstimv2(Pool,rate,raster,0,'twitter');
% 
%%
% close all
% CHANGE NAME %
% save('E:\DATA\ana_population\Dynamics_SUrate\PT_A1.mat','SUrate','-v7.3','-nocompression');
save('E:\DATA\combined\ana_population\Dynamics_SUrate\Trill_AL_3ani.mat','SUrate','-v7.3','-nocompression');
%% Extracting waveforms directly from units


N_list = [277,276,275,273,261,223];

for n = 1:length(N_list)
    rec_list = find([Pool{1}.neuron_nb] ==N_list(n));
    figure
    l_w = [0,0];
    
    for s = 1:2
        subplot(1,3,s)
        
        l_w(s) = size(Pool{1}(rec_list(s)).waveforms,2);
        list_w = randsample(1:l_w(s),100);
        for w = list_w
            plot(Pool{1}(rec_list(s)).waveforms(:,w))
            hold on
        end
        plot(mean(Pool{1}(rec_list(s)).waveforms,2),'Linewidth',3,'Color','k')
        hold off
        
        axis([0 60 -150 60])
    end
    X = [Pool{1}(rec_list(1)).waveforms Pool{1}(rec_list(2)).waveforms];
    W = pca(X);
    g = [ones(1,l_w(1)) ones(1,l_w(2))*2];
    g = g.';
    subplot(1,3,3)
    
    gscatter(W(:,1),W(:,2),g)
    drawnow
end













%% temporary code to add info to xbz files.
% 
% rmpath('Y:\xbz_Backup')
% 
% % addpath('\\datacenterchx.bme.jhu.edu\Project_MCR\xbz_Backup')
% addpath('D:\DATA\M12E\Experiments'); %path to xbz files
% 
% % i = 1;
% % M12Elog{1,2}
% 
% tic
% for i = 1:length(M12Elog)
%     if ~isempty(M12Elog{i,1})
%         if ~exist(['D:\DATA\M12E\Experiments\', M12Elog{i,1}, '.mat'])
%             x = eval(M12Elog{i,1});
%             x.hole_number = str2num(M12Elog{i,2});
%             x.track_number = str2num(M12Elog{i,3});
%             x.data = [];
%             save(['D:\DATA\M12E\Experiments\', M12Elog{i,1}, '.mat'],'x');
%         end
%     end
%     if mod(i,10) ==0
%         disp(M12Elog{i,1})
%         toc
%     end
% end
