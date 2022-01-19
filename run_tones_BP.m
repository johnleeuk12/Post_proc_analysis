

[Pool, rate, raster] =  Pool_BF([7,8,9,10],'M160E');

pause(0.1)
det_BF(Pool, rate, raster,'m160E',0);
%
out = ana_BF([7,8,9,10], 0,'m160E');


%% combine outs

out = {};

out.tags = out1.tags;
out.data{1} = [out1.data{1};out2.data{1}];

for on = 2:5
    out.data{on} = [out1.data{on};out2.data{on}];
end

holes = unique(out.data{1}(:,2));