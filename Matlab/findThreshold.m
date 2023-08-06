function [optiThres,F1] = findThreshold(predictor, csPCa)
% Function to find optimal threshold for defining cancer based on maximal
% F1-score
% positive = cancer
k=1;

for thres = min(predictor):0.1:max(predictor)
    xtemp = predictor>thres;
    P = sum(xtemp==1&csPCa==1)/sum(xtemp);
    R = sum(xtemp==1&csPCa==1)/sum(csPCa);
    F1(k) = (2*(P*R))/(P+R);
    threshold(k) = thres;
    k=k+1;
end

figure
plot(threshold, F1, '.k')
xlabel('threshold')
ylabel('F1')

[~ , idx] = max(F1);
optiThres = threshold(idx);

end