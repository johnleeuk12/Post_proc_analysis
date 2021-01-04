function [stim_label, stim_ticks] = stim_index(Pool,p)

code = Pool(p).xb.analysis_code;
ind = [];

switch code
    case 2320 %VT Phee CF
        ind = 10;
        ind2 = 1;
    case 2367 % Twitter CF
        ind = 9;
        ind2 = 1;
        
    case 2370 % Twitter IPI
        ind = 12;
        ind2 = 2;
        
    case 2335 % Trill CF
        ind = 10;
        ind2 = 1;
        
    case 2338 % Trill FM mode
    case 2341 %Trill AM mod
    case 2343 %Trill AM mod 
        ind = 33
        ind2 = 1;
    case 2339 %Trill FM rate
        ind = 22;
        ind2 = 1;
    
end

stim_label = Pool(p).xb.stimulus_ch1(:,ind);
for stim = 1:length(stim_label)
    %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
    stim_ticks{stim}=num2str(round(stim_label(stim)*10^ind2)/10^ind2);
end