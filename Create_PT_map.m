%% Code for tuning map

%{
The idea is to have a matrix that act as pixels on an image, with a
craniotomy being within a 7x7 grid. 
y axis negative is towards the LS
x axis positive is towards posterior

%}
h_positions = [1:15;-1*ones(1,15);-1*ones(1,15)].'; %positions of craniotomy center

%positions on x axis along the LS
h_positions(3,2) = 53;
h_positions(2,2) = h_positions(3,2) - 8;
h_positions(6,2) = h_positions(2,2) - 10;
h_positions(8,2) = h_positions(6,2) - 8;
h_positions(9,2) = h_positions(8,2) - 9;
h_positions(10,2) = h_positions(2,2)- 5;
h_positions(11,2) = h_positions(10,2) - 7;
h_positions(12,2) = h_positions(11,2) - 9;
h_positions(13,2) = h_positions(12,2) - 10;
h_positions(14,2) = h_positions(9,2) - 9;
h_positions(15,2) = h_positions(13,2) - 8;

%positions on y axis perpendicular to the LS
%position of LS = 1;
h_positions(3,3) = 3;
h_positions(2,3) = 4;
h_positions(6,3) = 5;
h_positions(8,3) = 5;
h_positions(9,3) = 5;
h_positions(10,3) = 10;
h_positions(11,3) = 13;
h_positions(12,3) = 10;
h_positions(13,3) = 11;
h_positions(14,3) = 4;
h_positions(15,3) = 10;

% positions of individual tracks.
t_positions = {}; %{x,y,bf,lat}




% hole 6
t_positions{6} = ...
    [[-2,2,12,18];... % relative track coordinates and best frequency
    [0,3,13,23];...
    [2,2,15,11];...
    [3,0,19,18];...
    ];


%hole 8
t_positions{8} = ...
    [[2,-1,10.5,18.4];...% relative track coordinates and best frequency
    [3,0,11.5,35];...
    [-1,-2,10.5,15.8];... 
    [0,0,8.5,27.7];...
    [-2,1,6.5,26.6];...
    ];


%hole 9
t_positions{9} = ...
    [[2,-1,4,26];...
    [0,0,5.4,14];...
    [-2,1,2.6,28];...
    [-2,-1,3.3,29];...
    
    ];

%hole 10
t_positions{10} = ...
    [[-2,2,7,22];... % relative track coordinates and best frequency
    [0,0,7,17];...
    [-2,-1,1.6,26];...
    [0,3,10,26];...
    [2,-1,4,47];...
    ];

%hole 11
t_positions{11} = ...
    [[0,-3,20.2,19.7];... 
    [3,0,13.3,24];...
    [-3,0,18.6,24];...
    [0,3,15.7,30];...
    [0,0,12.1,20];...
    ];

%hole 12
t_positions{12} = ...
    [[0,0,8.2,28.3];... 
    [-2,-1,10.4,15.3];...
    [0,3,11.8,29.1];...
    [-2,1,6.8,17.4];...
    ];

%hole 13
t_positions{13} = ...
    [[0,-3,4.2,23];... 
    [1,2,2.9,24.8];...
    [-3,0,5.4,18.8];...
    [3,0,7.4,30];...
    ];

%hole 15
t_positions{15} = ...
    [[0,0,2.3,36];... 
    [-2,-1,8.5,30];...
    [-1,-2,12.5,33];...
    ];


F = zeros(8*7,4*7);
L = zeros(8*7,4*7);
% plotting on to map F
for h = 1:15
    if ~isempty(t_positions{h})
        for t = 1:size(t_positions{h},1)
            coord = h_positions(h,2:3) + t_positions{h}(t,1:2);
            F(coord(1),coord(2)) = t_positions{h}(t,3);
            L(coord(1),coord(2)) = t_positions{h}(t,4);
        end
    end
end




F = F.';
% F = log2(F);    
% F_filt = imgaussfilt(F,0.5);

F_filt = log2(F);
% F_filt = F;
for s = 1:size(F_filt,1)*size(F_filt,2)
    if F_filt(s) == -Inf
        F_filt(s) = nan;
    end
end
figure
s= pcolor(F_filt);
s.EdgeColor = 'none';
caxis([0 5])

 cbh = colorbar ; %Create Colorbar
 cbh.Ticks = [0:5] ; %Create 8 ticks from zero to 1
 cbh.TickLabels = [{1},{2},{4},{8},{16},{32}];

 figure
 L = L.';
 for l = 1:size(L,1)*size(L,2)
     if L(l) ==0
         L(l) =nan;
     end
 end
 c = pcolor(L);
 c.EdgeColor = 'none';

 
% colormap(mycolormap)
 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    