function mu = calc_coherence(A)

% To calculate the maximum coherence value amongst the columns of sensing matrix A 
k = 1;
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if i < j
            in(k) = ((A(:,i))' * A(:,j)) / (sqrt(sum(abs(A(:,i)) .^ 2)) * sqrt(sum(abs(A(:,j)) .^ 2)));
            x(k) = i;
            y(k) = j;
            k = k + 1;
        end
    end
end
mu = max(abs(in));
end