function Create_tonotopy()

global D
D = [];
    animal_name = 'M160E';
    file_dir = 'D:\DATA\M160E\HT_position\';
    % Load EPhys Track positions
%%
    D.Ephys.FileCoorTracks = [file_dir, animal_name, ...
        '_track_positions.txt'];
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
    
    
    
    
    
%% Craniotomy positions

figure(1)

imshow(imread([file_dir, animal_name,'_map.png']));
f = imread([file_dir, animal_name,'_map.png']);
D.Crani.image_size.x = size(f,2);
D.Crani.image_size.y = size(f,1);


N_holes = length(D.Ephys.CoorTrackH);
for n = 1:N_holes
    h = drawpoint;
    D.Ephys.CoorTrackH{n}.hole_posi_x = h.Position(1);
    D.Ephys.CoorTrackH{n}.hole_posi_y = h.Position(2);
    pause(0.1)
    disp(['H',num2str(n),'... x position: ', num2str(h.Position(1)),...
        '      y position: ', num2str(h.Position(2))]);
    delete(h);
end
































