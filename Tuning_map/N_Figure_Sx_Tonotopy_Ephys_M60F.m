% function N_Figure_5_Tonotopy_Ephys
% 2 is for re-layout of panels (201227), delete TS full-view rotated panel.
% Figure Extended 6, Tonotopy, Through & Without Skull
%   Resolution analysis, R^2 on tuning & amp, across kernal sizes
%   Response temporal traces comparison
%   variance or SNR analysis (about repetition number)
%% Data Read in
global D G
D = [];     G = [];
% Gloabl Parameters
    D.Map_Ph =	300;
    D.Map_Pw =	480;
% Tonotopy, Stimuli Design
    D.AudTimeStim =     20;
    D.AudTimeStart =    2.7;
    D.AudTimeStop =     17.3;
    D.AudTimePip =      0.2;
    D.AudTimeSes =      400;
    % ColorMap
    addpath('D:\GitHub\fluffy-goggles\XinProc');
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
% Load Ephys Tuning Data
    D.Ephys.CF_H{4} = [];
    D.Ephys.H4_40 = load(    ['D:\Dropbox\==LightUp==\=Paper #4 Xintrinsic\Figure Code\'...
        'Data_Fig_Nx_M60F_H4_40dB.mat']);
    D.Ephys.H4_40 = D.Ephys.H4_40.X;
    for i = 1:length(D.Ephys.H4_40)
        h = D.Ephys.H4_40(i).Hole;
        t = D.Ephys.H4_40(i).Track;
        try
           D.Ephys.CF_H{h}.Tsu40{t} = [D.Ephys.CF_H{h}.Tsu40{t}  D.Ephys.H4_40(i).cf];
        catch
           D.Ephys.CF_H{h}.Tsu40{t} = D.Ephys.H4_40(i).cf;
        end
    end
    for h = 1:length(D.Ephys.CF_H)
        if isfield(D.Ephys.CF_H{h}, 'Tsu40')
            for t = 1:length(D.Ephys.CF_H{h}.Tsu40)
                D.Ephys.CF_H{h}.Tmean40(t) = mean(D.Ephys.CF_H{h}.Tsu40{t});
                D.Ephys.CF_H{h}.Tmean40st(t) = ...
                    log2(D.Ephys.CF_H{h}.Tmean40(t)/440)*12;
            end
        end
    end
% Load Ephys Tracks Data
    D.Ephys.FileCoorTracks = ['D:\Dropbox\==LightUp==\=Paper #4 Xintrinsic\Figure Code\'...
        'Data_Fig_Nx_M60F_map_tracks_xml_textcopied.txt'];
        % save the ppt with hole+track position slides only, as .xml.
        % reopen the .xml file in web browser(edge) 
        % resave the opened page into a .txt file
    D.Ephys.hFileCoorTracks = fopen(D.Ephys.FileCoorTracks);
    % read raw text
    i = 1;
    while ~feof(D.Ephys.hFileCoorTracks)
        D.Ephys.textline = fgetl(D.Ephys.hFileCoorTracks);
        D.Ephys.CoorTracksTextRaw{i} = D.Ephys.textline;      i = i+1;
    end
    fclose(D.Ephys.hFileCoorTracks);
    % extract useful text
    s=0;                Son=0;  % slide on;
    sp=0;   tSP = [];   SPon=0; % SP on;
    for i = 1:length(D.Ephys.CoorTracksTextRaw)
        D.Ephys.textline = D.Ephys.CoorTracksTextRaw{i};
        if Son ==0                      % search for a new slide
            if length(D.Ephys.textline)>=10
                if strcmp(D.Ephys.textline(1:10), '-<pkg:part') && ...
                   contains(D.Ephys.textline, 'officedocument.presentationml.slide+xml')
                    tt =    regexp(D.Ephys.textline, 'slide\d');
                    tt =	textscan(D.Ephys.textline(tt+5:end), '%d');
                    s = tt{1};  
                    Son = 1;    sp = 0;	% a new slide, with reset sp index
                end
            end
        else                            % work on the current slide
            if length(D.Ephys.textline)>=10
                if strcmp(D.Ephys.textline(1:10), '</pkg:part')	
                        Son = 0;        % close the slide
                end
            end
            if length(D.Ephys.textline)>=7
                if SPon == 0
                    if strcmp(D.Ephys.textline(1:7), '-<p:sp>')  
                    	SPon = 1;   tSP = [];	% start a new "sp"
                    end
                else
                    if strcmp(D.Ephys.textline(1:7), '</p:sp>')  
                        sp = sp + 1;            % finish and save the current "sp"
                        D.Ephys.CoorTracksS{s}.SP{sp} = tSP;
                    	SPon = 0;               % close the current "sp"
                    end
                    if length(D.Ephys.textline)>=10 % fill up the "sp" information
                        if      strcmp(D.Ephys.textline(2:9), '<p:cNvPr')% p:cNvPr Name
                            tSP.type = D.Ephys.textline(17:20);
                        elseif	strcmp(D.Ephys.textline(1:5), '<a:t>')	% Text
                            tSP.text = D.Ephys.textline(6:end);
                        elseif  strcmp(D.Ephys.textline(1:6), '<a:off')	% off
                            tSP.off =  D.Ephys.textline(8:end);
                        elseif  strcmp(D.Ephys.textline(1:6), '<a:ext')	% ext
                            tSP.ext =  D.Ephys.textline(8:end);
                        end
                    end
                end
            end
        end
    end
    % Construct the tracks of each hole
    for s = 1:length(D.Ephys.CoorTracksS)
        for sp = 1:length(D.Ephys.CoorTracksS{s}.SP)
            switch D.Ephys.CoorTracksS{s}.SP{sp}.type
                case 'Titl'     % Hole #
                    h = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.text(2:end), '%d');      h = h{1};
                        D.Ephys.CoorTrackH{h} = [];   
                case 'Oval'     % Hole Contour
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.off, 'x="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.off(tt+3:end),'%d');    tt = tt{1};
                        D.Ephys.CoorTrackH{h}.holeoffx = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.off, 'y="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.off(tt+3:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.holeoffy = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.ext, 'cx="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.ext(tt+4:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.holeextx = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.ext, 'cy="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.ext(tt+4:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.holeexty = double(tt);                    
                case 'Text'     % Each Track in the hole
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.text(2:end),'%d');      t = tt{1};
                        D.Ephys.CoorTrackH{h}.track{t} = [];    
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.off, 'x="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.off(tt+3:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.track{t}.offx = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.off, 'y="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.off(tt+3:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.track{t}.offy = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.ext, 'cx="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.ext(tt+4:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.track{t}.extx = double(tt);
                    tt = regexp(D.Ephys.CoorTracksS{s}.SP{sp}.ext, 'cy="\d');
                    tt = textscan(D.Ephys.CoorTracksS{s}.SP{sp}.ext(tt+4:end),'%d');	tt = tt{1};
                        D.Ephys.CoorTrackH{h}.track{t}.exty = double(tt);                          
            end
        end
    end
% Load Ephys Holes Data
    % Surface Picture:  <p:cNvPr id="4" name="Picture 3">
    %                   <a:off x="609599" y="0"/>
    %                   <a:ext cx="10972801, 'exty', 6858000"/>
    D.Ephys.CoorPic.offx = 609599;
    D.Ephys.CoorPic.offy = 0;
    D.Ephys.CoorPic.extx = 10972801;
    D.Ephys.CoorPic.exty = 6858000;
    D.Ephys.CoorPic.angle = 18;    
            % the background pic is 18 degree clockwise than the ephys pic
            % or the ephys pic is positioned as 342 degree (in ppt)
        % <a:off x="2559904" y="3316033"/>    <a:ext cx="510174" cy="485506"/>    <a:t>2</a:t>
        % <a:off x="3090019" y="3163562"/>    <a:ext cx="510174" cy="485506"/>    <a:t>1</a:t>
        % <a:off x="2868699" y="3672318"/>    <a:ext cx="510174" cy="485506"/>    <a:t>3</a:t>
        % <a:off x="3209508" y="4106247"/>    <a:ext cx="510174" cy="485506"/>    <a:t>4</a:t>
        % <a:off x="3465571" y="2976686"/>    <a:ext cx="510174" cy="485506"/>    <a:t>5</a:t>
        % <a:off x="3841123" y="2789810"/>    <a:ext cx="510174" cy="485506"/>    <a:t>6</a:t>
        % <a:off x="2369403" y="3775424"/>    <a:ext cx="510174" cy="485506"/>    <a:t>7</a:t>
        % <a:off x="2739528" y="4116270"/>    <a:ext cx="510174" cy="485506"/>    <a:t>9</a:t>
        % <a:off x="2559904" y="4313830"/>    <a:ext cx="510174" cy="485506"/>    <a:t>8</a:t>
        % <a:off x="2972585" y="4533909"/>    <a:ext cx="510174" cy="485506"/>    <a:t>10</a:t>
        % <a:off x="3424419" y="4575696"/>    <a:ext cx="510174" cy="485506"/>    <a:t>11</a:t>
    D.Ephys.CoorHole{1} =   struct('offx', 3090019, 'offy', 3163562, 'extx', 510174, 'exty', 485506); %1</a:t>
    D.Ephys.CoorHole{2} =   struct('offx', 2559904, 'offy', 3316033, 'extx', 510174, 'exty', 485506); %2</a:t>
    D.Ephys.CoorHole{3} =   struct('offx', 2868699, 'offy', 3672318, 'extx', 510174, 'exty', 485506); %3</a:t>
    D.Ephys.CoorHole{4} =   struct('offx', 3209508, 'offy', 4106247, 'extx', 510174, 'exty', 485506); %4</a:t>
    D.Ephys.CoorHole{5} =   struct('offx', 3465571, 'offy', 2976686, 'extx', 510174, 'exty', 485506); %5</a:t>
    D.Ephys.CoorHole{6} =   struct('offx', 3841123, 'offy', 2789810, 'extx', 510174, 'exty', 485506); %6</a:t>
    D.Ephys.CoorHole{7} =   struct('offx', 2369403, 'offy', 3775424, 'extx', 510174, 'exty', 485506); %7</a:t>
    D.Ephys.CoorHole{8} =   struct('offx', 2559904, 'offy', 4313830, 'extx', 510174, 'exty', 485506); %8</a:t>
    D.Ephys.CoorHole{9} =   struct('offx', 2739528, 'offy', 4116270, 'extx', 510174, 'exty', 485506); %9</a:t>
    D.Ephys.CoorHole{10} =  struct('offx', 2972585, 'offy', 4533909, 'extx', 510174, 'exty', 485506); %10</a:t>
    D.Ephys.CoorHole{11} =  struct('offx', 3424419, 'offy', 4575696, 'extx', 510174, 'exty', 485506); %11</a:t>
% calculate the relative track position to the hole center
    for h = 1:length(D.Ephys.CoorTrackH)
        D.Ephys.CoorHole{h}.centerxy = [...
            D.Ephys.CoorHole{h}.offx + D.Ephys.CoorHole{h}.extx/2,...
            D.Ephys.CoorHole{h}.offy + D.Ephys.CoorHole{h}.exty/2];
        D.Ephys.CoorHole{h}.radiusxy = [...
            D.Ephys.CoorHole{h}.extx   D.Ephys.CoorHole{h}.exty]/2;
        if ~isempty(D.Ephys.CoorTrackH{h})
            hcenterxy = [   D.Ephys.CoorTrackH{h}.holeoffx + D.Ephys.CoorTrackH{h}.holeextx/2 ...    
                            D.Ephys.CoorTrackH{h}.holeoffy + D.Ephys.CoorTrackH{h}.holeexty/2   ];
            for t = 1:length(D.Ephys.CoorTrackH{h}.track)
                trackxy = [ D.Ephys.CoorTrackH{h}.track{t}.offx + D.Ephys.CoorTrackH{h}.track{t}.extx/2 ...
                            D.Ephys.CoorTrackH{h}.track{t}.offy + D.Ephys.CoorTrackH{h}.track{t}.exty/2 ]; 
                trackxyrela = (trackxy - hcenterxy)./...
                            ([D.Ephys.CoorTrackH{h}.holeextx D.Ephys.CoorTrackH{h}.holeexty]/2);
                [tracktrrela(1), tracktrrela(2)] = cart2pol(trackxyrela(1), trackxyrela(2)); 
                                                                    % theta angle is clockwise here due to the reverse y
                tracktrrelaccw = [-tracktrrela(1) tracktrrela(2)];  % theta angle counter-clockwise
                D.Ephys.CoorTrackH{h}.track{t}.trackxy =        trackxy;
                D.Ephys.CoorTrackH{h}.track{t}.trackxyrela =    trackxyrela;
                D.Ephys.CoorTrackH{h}.track{t}.tracktrrela =    tracktrrela;
                D.Ephys.CoorTrackH{h}.track{t}.tracktrrelaccw = tracktrrelaccw;
            end
        end
    end
% plot     
    figure;
    % the map background
    plot(   D.Ephys.CoorPic.offx + D.Ephys.CoorPic.extx * [0 1 1 0 0],...
            D.Ephys.CoorPic.offy + D.Ephys.CoorPic.exty * [0 0 1 1 0]);
    hca = gca;
    hca.NextPlot =  'add';
    hca.YDir =      'reverse';
    hca.Color =     0*[1 1 1];
    hca.XLim =      D.Ephys.CoorPic.offx + D.Ephys.CoorPic.extx *[-0.1 1.1];
    hca.YLim =      D.Ephys.CoorPic.offy + D.Ephys.CoorPic.exty *[-0.1 1.1];
    % holes & tracks
    for h = 1:length(D.Ephys.CoorTrackH)
        % holes
        plot(   D.Ephys.CoorHole{h}.offx + D.Ephys.CoorHole{h}.extx*(cos(linspace(0,2*pi, 181))+1)/2,...
                D.Ephys.CoorHole{h}.offy + D.Ephys.CoorHole{h}.exty*(sin(linspace(0,2*pi, 181))+1)/2,...
                    'Color',                0.7*[1 1 1],...
                    'LineWidth',            0.75);
        if ~isempty(D.Ephys.CoorTrackH{h})
            for t = 1:length(D.Ephys.CoorTrackH{h}.track)
                [pictrackxyrela(1), pictrackxyrela(2)] = pol2cart(...
                    D.Ephys.CoorTrackH{h}.track{t}.tracktrrela(1) - D.Ephys.CoorPic.angle/180*pi,...
                  	D.Ephys.CoorTrackH{h}.track{t}.tracktrrela(2) );
                pictrackxy = D.Ephys.CoorHole{h}.centerxy + ...
                    pictrackxyrela.*D.Ephys.CoorHole{h}.radiusxy;
                % tracks
                text( pictrackxy(1), pictrackxy(2), sprintf('T%d', t),...
                    'Color',                0.25*[1 1 1],...
                    'Rotation',             D.Ephys.CoorPic.angle,...
                    'HorizontalAlignment',	'Center',...
                	'VerticalAlignment',	'Middle');
                try
%                     scatter(pictrackxy(1), pictrackxy(2), 36,...
%                         D.Ephys.CF_H{h}.Tmean40st(t), 'filled',...
%                         'MarkerFaceAlpha',       0.7);
                    scatter(pictrackxy(1), pictrackxy(2), 900,...
                        D.Ephys.CF_H{h}.Tmean40st(t), 'filled',...
                        'MarkerFaceAlpha',       0.7);
                end
            end
        end
    end
%         colormap('jet');
                colormap(squeeze(G.A_CM_RGB(1,:,:)))
                caxis([0 72])
    
    return;