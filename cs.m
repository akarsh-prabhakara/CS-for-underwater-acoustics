function Xest = cs(theta,M,ratio,L,snr,W)

% Initialisation
deg = -90:1:90;

% Initialising X according 
x = (randn(length(theta),L)) .* (sqrt(0.5)) + 1j * (randn(length(theta),L)) .* (sqrt(0.5));
X = zeros(length(deg),L);
for i = 1:length(deg)
    for j = 1:length(theta)
        if deg(i) == theta(j)
            X(i,:) = x(j,:);
        end
    end
end

% Forming the steering vector for each source
sensorN = 0:M-1;
%sensorN = 0 + (M-1)*rand(M,1);
for i = 1:length(deg)
    A(:,i) = exp(1j * 2 * pi * sind(deg(i)) * sensorN * ratio)';
end

% Obtaining the noise free sensor measurements for L snapshots
for i = 1:L
    AX(:,i) = A * X(:,i);
end

% Adding AWGN
for i = 1:L
    n(i) = (sum(AX(:,i) .* conj(AX(:,i))))^0.5;
    N(:,i) = ((randn(length(sensorN),1))*(sqrt(0.5*((n(i)*10^(-snr/20))^2)/M))) +1j * ((randn(length(sensorN),1))*(sqrt(0.5*((n(i)*10^(-snr/20))^2)/M)));
end     
Y = (AX + N);

for i = 1:size(Y,1)
    Yl2(i,1) = mean(Y(i,:));
end
%for i =1:L   
cvx_begin quiet
    variable xest(length(deg),1)
    minimize(norm(diag(W)*xest,1)) 
    subject to
        norm(Yl2-A*xest,2)<=0.8*mean(n)
cvx_end
Xest = xest;
%end
end