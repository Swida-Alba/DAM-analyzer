function [varargout] = cosinor(Y,t,P)
    arguments
        Y; % activity data, 1 x N vector
        t; % in hour % time data corresponding to Y, each element at i is the time of the corresponding y_i in Y
        P = 24; % in hour % given period
    end
    N = length(Y);
    Y_bar = mean(Y);
    x = cos(2*pi*t/P);
    z = sin(2*pi*t/P);
    d = [sum(Y);
        sum(Y.*x);
        sum(Y.*z)];
    S = [    N,    sum(x),    sum(z);
        sum(x), sum(x.^2), sum(x.*z);
        sum(z), sum(x.*z), sum(z.^2)];
    u_est = S^(-1) * d; % u estimate
    
    M_est     = u_est(1); % MESOR, Midline Estimating Statistic of Rhythm
    beta_est  = u_est(2);
    gamma_est = u_est(3);

    A = (beta_est ^ 2 + gamma_est ^ 2) ^ (1/2); % amplitude estimate
    phi = mod(atan(-gamma_est / beta_est),pi); % acrophase estimate
    Y_est = M_est + beta_est * x + gamma_est * z; % estimated activity
    RSS = sum((Y - Y_est) .^ 2);
    MSS = sum((Y_est - Y_bar) .^ 2);
    F = (MSS/2) / (RSS/(N-3));

    p = 1 - fcdf(F,2,N-3); % p-value of F-test
    stat = struct('p', p, 'F', F, 'RSS', RSS, 'MSS', MSS);

    varargout = {stat, M_est, A, phi, beta_est, gamma_est};

end