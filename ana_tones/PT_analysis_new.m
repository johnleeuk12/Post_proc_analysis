

nH = 15;
nTr = 2;


C.H{1,nH}{1,nTr}.BF  
C3.H{1,nH}{1,nTr}.BF  

O = [];
p = 1;
for nH = 1:16
    if ~isempty(C.H{1,nH})
        for nTr = 1:size(C.H{1,nH},2)
            if ~isempty(C.H{1,nH}{1,nTr})
                O(1,p) = C.H{1,nH}{1,nTr}.BF.mean;
                O(2,p) = C3.H{1,nH}{1,nTr}.BF.mean;
            end
            p = p+1;
        end
        
    end
end
          

scatter(O(1,:),O(2,:))