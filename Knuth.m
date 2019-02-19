function L = Knuth(lambda)
    L = exp(-lambda);
    k = 0;
    p = 1;
    while p > L
        k = k + 1;
        p = p*rand(1);
    end
    L = k - 1;