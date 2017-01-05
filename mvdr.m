function [P,ang] = mvdr(theta,M,ratio,L,snr)

% Initialisation
res = 181;
deg = linspace(-90,90,res);

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
Y = AX + N;

% Forming the steering vectors for 1D search

sensorN = 0:M-1;
for k = 1:length(deg)
    Z(:,k) = exp(1j * 2 * pi * sind(deg(k)) * sensorN * ratio)';
end

% Estimating the auto correlation matrix
R = zeros(size(Y,1),size(Y,1));
for i = 1:L
    R = R + Y(:,i) * (Y(:,i))';
end
R = R / L;

for k = 1:length(deg)
    W = ((R) \ Z(:,k)) / ((Z(:,k))' * ((R) \ Z(:,k)));
    P(k) = (W' * R * W);
end

[P_sort,loc] = findpeaks((abs(P)));
[~,pos] = sort(P_sort);
ang = sort(deg(loc(pos(length(P_sort)-length(theta)+1:end))));
end