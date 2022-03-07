

[Pool, rate, raster] =  Pool_BF([1:5],'M56E');

pause(0.1)
det_BF(Pool, rate, raster,'m56E',0);
%
out = ana_BF([1:5], 0,'m56E');


%% combine outs

out = {};

out.tags = out1.tags;
out.data{1} = [out1.data{1};out2.data{1}];

for on = 2:5
    out.data{on} = [out1.data{on};out2.data{on}];
end

holes = unique(out.data{1}(:,2));

%% noise combine noise 

out2 = ana_noise([7:10],'M160E',B,N_list);