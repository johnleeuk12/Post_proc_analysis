%% using comp from ana_singleunits for additional analysis and figure generation 04/2022
div_ind = 4;

peak_ind4 = comp.peak_ind3*1/div_ind;
DR_scatter = comp.DR_peak_scatter;



%% stats for comp


M = zeros(1,20);
for d = 2:21
    ind = [find(comp.peak_ind3 == d-1); find(comp.peak_ind3 == d); find(comp.peak_ind3 == d+1)];
%     ind = find(peak_ind3 ==d);
    M(d) = median(comp.DR_peak_scatter(ind,:));
end


figure(59)
hold off

scatter(comp.peak_ind3*1/div_ind, comp.DR_peak_scatter);
axis([0,5,-0.05,1.05])
hold on
plot([0:1/div_ind:5],M,'LineWidth',2);
plot([0:1/div_ind:5],fit([0:1/div_ind:5]),'--g','LineWidth',2);

figure
for a = [ -.4947 -.5781]
% a = -0.5;
x = [0:0.25:5];
eqn = @(x) -1*exp(a*x)+1;
plot(x,eqn(x),'Linewidth',2)
hold on

end

% std_scatter = zeros(1,10);
% for p_ind = 1:10
%     s_ind = find(ceil(comp.peak_ind3/2) == p_ind);
%     std_scatter(p_ind) = std(comp.DR_peak_scatter(s_ind,1));
% end
% 
% plot([1:10]/2,std_scatter)
%% combining comp files
%% combining comp 
comp1 = load('E:\DATA\Combined\ana_single_units\Trill\Trill_AL_HF_0SD.mat');
comp2 = load('E:\DATA\Combined\ana_single_units\Trill\Trill_AL_LF_0SD.mat');

comp = {};
comp.DR_peak = [comp1.comp.DR_peak;comp2.comp.DR_peak];
comp.g= [comp1.comp.g;comp2.comp.g];
comp.peak_ind3 = [comp1.comp.peak_ind3;comp2.comp.peak_ind3];
comp.DR_peak_scatter = [comp1.comp.DR_peak_scatter;comp2.comp.DR_peak_scatter];

div_ind = 2;

peak_ind4 = comp.peak_ind3*1/div_ind;
DR_scatter = comp.DR_peak_scatter;


%% comparing comp files 3/2022

group1 = comp1.comp.DR_peak(find(comp1.comp.g == 2),2);
group2 = comp2.comp.DR_peak(find(comp2.comp.g == 2),2);

p = ranksum(group1,group2)

figure(4)
edges = 0:0.05:1;
histogram(group1,edges);
hold on
histogram(group2,edges);









%% testing overlapping units
overlap_ind = false(length(comp.Sel_nid_all),1);
overlap_ind1 = false(length(comp1.Sel_nid_all),1);


for nn = 1:length(comp1.Sel_nid_all)
    ind = find(comp.Sel_nid_all == comp1.Sel_nid_all(nn) & comp.Sel_aid_all == comp1.Sel_aid_all(nn));
    if ~isempty(ind)
        overlap_ind(ind) = true;
        overlap_ind1(nn) = true;
    end
end

nan_ind1 = ~isnan(comp.g);
overlap_ind = overlap_ind(nan_ind1);
figure(61)
scatter(comp.peak_ind3*1/div_ind, comp.DR_peak_scatter)
hold on
scatter(comp.peak_ind3(overlap_ind)/div_ind,comp.DR_peak_scatter(overlap_ind))
axis([0,5,-0.05,1.05])
