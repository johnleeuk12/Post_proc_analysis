for vc = 1:2:9
    [DD, good_list_nid,good_list,aid_list] = format4DataHigh(SUrate,10,1);
    
    PreStim = 300;
    PostStim = 500;
    nreps = 10;
    
    
    
    VT_code = 22;
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
    
    ana_corr2
        drawnow
    disp(vc)
    pause()

end