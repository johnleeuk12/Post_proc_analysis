function create_map()
%% map out analysis results using C and D
animal_name = 'M60F';

%rotation
% if right hemisphere
theta = -27;
% if left hemisphere
%         theta = 27;

% addpath('D:\GitHub\Post_proc_analysis\util');
load(fullfile('E:\DATA', filesep, animal_name, filesep,'HT_position\D_tonotopy.mat'));
file_dir = fullfile('E:\DATA', filesep, animal_name, filesep,'HT_position\');

global D C

%% Create track positions for stim panels 03/26/2022
% comment out when plotting tonotopy figures

% mat = [
% 1	2
% 1	4
% 1	5
% 2	5
% 2	8
% 3	2
% 3	3
% 3	5
% 3	7
% 5	2
% 5	3
% 5	5
% 
% ];
% 
% h_list = unique(mat(:,1));
% 
% C = {};
% for h = 1:length(h_list)
%     t_list = unique(mat(find(mat(:,1) == h_list(h)),2));
%     for tr = 1:length(t_list)
%         %         freq_range = mat(find(mat(:,1) == h_list(h)),3);
%         
%         %     C.H{1,h_list(h)}{1,t_list(tr)} = freq_range(tr);
%         C.H{1,h_list(h)}{1,t_list(tr)} = 1;
%     end
% end
%     

%% colormap

G.A_CM_TuneHueTplt =     HueRedist('HSVadjusted', 'Circular');
G.A_CM_TuneHueTpltIn =   linspace(0, 1, length(G.A_CM_TuneHueTplt))';
G.A_CM_AmpNum =     100;
G.A_CM_GrayNum =	0;
G.A_CM_HSV(:,:,1) =	repmat(G.A_CM_TuneHueTplt', G.A_CM_AmpNum+1, 1);
G.A_CM_HSV(:,:,2) =	repmat(linspace(1, 0, G.A_CM_AmpNum+1)', 1, length(G.A_CM_TuneHueTplt));
G.A_CM_HSV(:,:,3) =	repmat(linspace(1, 0, G.A_CM_AmpNum+1)', 1, length(G.A_CM_TuneHueTplt));
G.A_CM_RGB =  reshape(...
    hsv2rgb( reshape(G.A_CM_HSV,[],3) ),...
    G.A_CM_AmpNum+1, length(G.A_CM_TuneHueTplt), 3);



%%
N_holes = length(D.Ephys.CoorTrackH);
    

    
figure(1)
hold off


f = imread([file_dir, animal_name,'_chb2_with_C.png']);
imshow(f);
% imshow(ones(size(f,1),size(f,2),3))
hold on
% creating a tonotopy map

for n = 1:N_holes
    center = zeros(1,2);
    try
        center = [D.Ephys.CoorTrackH{n}.holeoffx + D.Ephys.CoorTrackH{n}.holeextx/2, ...
                    D.Ephys.CoorTrackH{n}.holeoffy + D.Ephys.CoorTrackH{n}.holeexty/2];
        trackposi = zeros(length(D.Ephys.CoorTrackH{n}.track),2);
        
        for tr = 1:length(D.Ephys.CoorTrackH{n}.track)
            trackposi(tr,:) = [D.Ephys.CoorTrackH{n}.track{1,tr}.offx + D.Ephys.CoorTrackH{n}.track{1,tr}.extx/2, ...
                                 D.Ephys.CoorTrackH{n}.track{1,tr}.offy + D.Ephys.CoorTrackH{n}.track{1,tr}.exty/2];
        end
        posi = (trackposi-center)/D.Ephys.CoorTrackH{n}.holeextx;
%         posi(:,2) = -posi(:,2);

        R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
        posi = posi*R';
        
        W = 50;
        posi = posi*W + [D.Ephys.CoorTrackH{n}.hole_posi_x, D.Ephys.CoorTrackH{n}.hole_posi_y]; %width of craniotomy
        
        % plotting wanted variable
        if ~isempty(C.H{1,n})
            for tr = 1:length(C.H{1,n})
                if ~isempty(C.H{1,n}{1,tr})
                    % use this to check number of tracks for a given stim
                    V = C.H{1,n}{1,tr}.BF.mean;
%                     scatter(posi(tr,1),posi(tr,2),100,V,'k');
                    % end
                    
%                     V = C.H{1,n}{1,tr}.minlat.med-200; %-C.H{1,n}{1,tr}.minlat.med;%/(32*1e3); % change here to change variable
                    
                                        V = log2(V/440)*12;
%                     V = log2(V*1e-3);
                    
                    scatter(posi(tr,1),posi(tr,2),100,V,'filled','MarkerFaceAlpha',0.7);
                    
                end
                
            end
        end
        
        
    catch
        fprintf(['H %4d does not have tracks \n'],n')
    end
end
% 
colormap(flipud(parula))
colormap(squeeze(G.A_CM_RGB(1,:,:)))

% colormap(parula)
% caxis([0 1]);
caxis([0 72]);
colorbar




% caxis([0 1]);
% caxis([-1 5])

