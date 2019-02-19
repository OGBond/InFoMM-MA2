N = 3e3;            % Number of observations
K = 3;              % Number of clusters to be used

param = randi(100,1,K);  % Randomly generated Poisson parameters 
data = zeros(1,N);

% Generates N "pseudo-data" where each parameter to use for a Poisson
% distribution is randomly selected to generate an outcome
for i = 1:N
    data(i) = poissrnd(param(randi(K)));
end

tau = zeros(N,K);   % Storage space for latent variables
F = zeros(N,K);     % Storage space for probability distributions

lambda = max(param)*rand(1,K);   % Storage space for initial random guess of lambdas
pi = rand(1,K);                  % Storage space for initial random guess of mixture coefficients

% Generate initial random distribution of latent variables
for i = 1:N
    probs = rand(1,K);
    tau(i,:) = probs/norm(probs,1);
end

for j = 1:1e3
    
    % Updates values of pi and lambda guesses according to maximum
    % likelihood estimatiors
    
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
    
    % Updates matrix of soft latent variables
    
    for i = 1:N
        for k = 1:K
            tau(i,k) = (pi(k).*F(i,k)) ./ (F(i,:)*pi');
        end
    end
end

histogram(data,'FaceColor','y','Normalization','probability')
hold on

P1 = poisspdf(0:max(data),lambda(1));
P2 = poisspdf(0:max(data),lambda(2));
P3 = poisspdf(0:max(data),lambda(3));

plot(P1,'r')
plot(P2,'g')
plot(P3,'m')

legend(["Original data";
        "\lambda = " + num2str(lambda(1));
        "\lambda = " + num2str(lambda(2));
        "\lambda = " + num2str(lambda(3))])

title('Visualisation of data consisting of three Poisson distributions')
xlabel('Data values')
ylabel('Probability density')