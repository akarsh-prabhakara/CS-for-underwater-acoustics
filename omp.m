function [s,sig,loc] = omp(A,y)
r = {y};
set = [];
k = 1;
for t = 1:3
%while sum(abs(r{k}).^2) > sum(abs(y).^2)*10^-7
    
    ip = (A' * r{k}).' ./ (sqrt(sum(abs((A)) .^ 2)));
    %ip = (A' * r{k}) ./ sum(A .^ 2);
    [~,loc(k)] = max(abs(ip));
    set(:,k) = A(:,loc(k));
    P = set * inv(set' * set) * set';
    r{k+1} = (eye(size(P,1)) - P) * y;
    k = k + 1;
end
sig = pinv(set) * (P * y);
s = zeros((size(A,2)),1);
s(loc) = sig;
end
