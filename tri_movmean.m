function O = tri_movmean(X2,S)
%%
% function for moving window mean with octave range
% S is window size, and moving window is S/2
% X must be an array of 2 columns, with the first column CF and second
% anything.
% S = 0.125;

range = [2,32];

fstart = range(1);

edges = fstart;
while fstart <= range(2)
    fend = fstart*2^S;
    edges = [edges fend];
    fstart = fend;
end
    
% calculating average in a moving window
B = X2(:,1);

O = zeros((length(edges)-2),2);
for i = 1:length(edges)-2
    ind = find(B<edges(i+2) & B>edges(i));
    weight = zeros(1,length(ind));
    for f = 1:length(ind)
        if B(ind(f))<edges(i+1)
            weight(f) = 1- log2(edges(i+1)/B(ind(f)))/S;
        elseif B(ind(f)) >=edges(i+1)
            weight(f) = 1- log2(B(ind(f))/edges(i+1))/S;        
        end
    end
    O(i,2) = mean(X2(find(B<edges(i+2) & B>edges(i)),2).*weight.');
    O(i,1) = edges(i+1);
end
    





















end