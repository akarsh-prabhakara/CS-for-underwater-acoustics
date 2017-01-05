function [Xest] = cs_rw1(theta,M,ratio,L,snr,sigma)
deg = -90:1:90;
sensorN = 0:M-1;
x = (randn(length(theta),L))*sqrt(0.5)+1j*(randn(length(theta),L))*sqrt(0.5);
X = zeros(length(deg),L);
%for i = 1:length(deg)
%    for j = 1:length(theta)
%        if deg(i) == theta(j)
%            X(i,:) = x(j,:);
%        end
%    end
%end
for i = 1:length(theta)
    A(:,i) = exp(1j * 2 * pi * sind(theta(i)) * sensorN * ratio)';
end
for i = 1:L
    Ax(:,i) = A * x(:,i);
end
% sensorN = 0 + (M-1)*rand(M,1);
for i = 1:L
    n(i) = (sum(Ax(:,i) .* conj(Ax(:,i))))^0.5;
    N(:,i) = randn(M,1)*sqrt(0.5*(n(i)*10^(-snr/20))^2/M)+1j*randn(M,1)*sqrt(0.5*(n(i)*10^(-snr/20))^2/M);
    no(i) = norm(N(:,i));
end    
Y = Ax + N;
%for i = 1:length(deg)
%    A(:,i) = exp(1j * 2 * pi * sind(deg(i)) * sensorN * ratio)';
%end
for i = 1:length(deg)
    A(:,i) = exp(1j * 2 * pi * sind(deg(i)) * sensorN * ratio)';
end

if L == 1
i = 0;
xest = ones(length(deg),1);
old_est = ones(length(deg),1)*20;
W = ones(length(deg),1);
%sigma = 0.0001;
while (old_est-xest).^2 >= 10^-8
i = i + 1;
if i == 1
cvx_begin quiet
    variable xest(length(deg),1) complex
    minimize(norm(diag(W)*xest,1)) 
    subject to
        norm(Y-A*xest,2)<=norm(N)
cvx_end
else
old_est = xest;    
for j = 1:length(deg)
W(j) = 1 / (abs(old_est(j))+sigma);
end
cvx_begin quiet
    variable xest(length(deg),1) complex
    minimize(norm(diag(W)*xest,1)) 
    subject to
        norm(Y-A*xest,2)<=norm(N)
cvx_end
end
%plot(abs(xest))
%pause(0.1)
end
else
i = 0;
Xcap = ones(length(deg),1);
old_est = ones(length(deg),1)*20;
W = ones(length(deg),1);
%sigma = 10;
while norm(old_est-Xcap,'fro') >= 10^-8
i = i + 1;
if i == 1
cvx_begin
    variable Xcap(size(X,1),size(X,2)) complex
    minimize(sum(W.*(sum((Xcap.*conj(Xcap)),2)))+norm(Y-A*Xcap,'fro')) 
    %subject to
    %    norm(Y-A*Xcap,'fro')<=mean(no)  
cvx_end
else
old_est = Xcap;    
for j = 1:length(deg)
W(j) = 1 / (abs(old_est(j))+sigma);
end
cvx_begin
    variable Xcap(size(X,1),size(X,2)) complex
    minimize(sum(W.*(sum((Xcap.*conj(Xcap)),2)))+norm(Y-A*Xcap,'fro')) 
    %subject to
    %    norm(Y-A*Xcap,'fro')<=mean(no) 
cvx_end
end    
end
end
Xest = xest;
end
