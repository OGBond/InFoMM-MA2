data =  load("PoissonData.txt");    % Loading Pete Grindrod's data set

N = length(data);   % Number of observations
K = 2;              % Number of clusters

tau = zeros(N,K);   % Storage space for latent variables
F = zeros(N,K);     % Storage space for probability distributions

lambda = 10*rand(1,K);   % Storage space for Possion parameters
pi = rand(1,K);       % Storage space for mixture coefficients

% Generate initial random distribution of latent variables
for i = 1:N
    probs = rand(1,K);
    tau(i,:) = probs/norm(probs,1);
end

while abs(min(max(tau)) - 1) > 1e-9
    
    % Updates values of pi and lambda
    
    for k = 1:K   
        sum = 0;
        for i = 1:N
            sum = sum + data(i)*tau(i,k);
        end
        pi(k) = norm(tau(:,k),1)/N;
        lambda(k) = sum/(N*pi(k));
    end
    
    % Creates values of log(F) which is a Poisson distribution
    % mass function
    
    for i = 1:N
        for k = 1:K
            logF = -lambda(k) + data(i)*log(lambda(k)) - log(factorial(data(i)));
            F(i,k) = exp(logF);
        end
    end
    
    for i = 1:N
        for k = 1:K
            tau(i,k) = (pi(k).*F(i,k)) ./ (F(i,:)*pi');
        end
    end
end

histogram(data,'FaceColor','y','Normalization','probability')
hold on

P1 = poisspdf(0:80,lambda(1));
P2 = poisspdf(0:80,lambda(2));

plot(P1,'r')
plot(P2,'g')

legend(["Original data";
        "\lambda = " + num2str(lambda(1));
        "\lambda = " + num2str(lambda(2))])
    
xlim([0,80])
title('Visualisation of data consisting of two Poisson distributions')
xlabel('Data values')
ylabel('Probability density')