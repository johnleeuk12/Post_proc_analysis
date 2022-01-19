
threshX3 = mean2(X(:,1:150))+std2(X(:,1:150))*2;
min_lat_pool = [];
if ~isempty(BF_ind)
    if BF_ind < 3
        X3 = X(BF_ind:BF_ind+2,200:end);
    elseif BF_ind > size(X,2)-3
        X3 = X(BF_ind-2:BF_ind,200:end);
    else
        X3 = X(BF_ind-2:BF_ind+2,200:end);
    end
    
%     thresholding method 1
    
    X3 = mean(X3,1); % thresholding method 2 Recanzone et al 2000. comment to not use
    
    
    for nx = 1:size(X3,1)
        t = 1;
        tt = 1; %counter to break out of while loop
        while (mean(X3(nx,t:t+1))< threshX3 || mean(X3(nx,t+2:t+3))< threshX3 || mean(X3(nx,t+4:t+5))< threshX3) && tt~=0
            t = t+1;
            if t > size(X3,2)-6
                pause
                disp(t)
                pause
                tt = 0;
                break
            end
        end
        if tt ~=0
            min_lat_pool = [min_lat_pool, t+200];
        end
    end
    
end
