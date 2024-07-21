function [varargout] = CircadianPeriodogram(Activity, binWidth, PeriodRange, options)
% Enright (chi-square) is from Sokolove P.G., & Bushell W.N., 1978, J Theor Biol; Refinetti R., ..., & Halberg G., 2007, Biol Rhythm Res.
% Lomb-Scargle is from Refinetti R., ..., & Halberg G., 2007, Biol Rhythm Res.
% Cosinor is from Cornelissen G., 2014, Theor Biol Med Model; Refinetti R., ..., & Halberg G., 2007, Biol Rhythm Res;
% greedy chi-square is from Tackenberg M.C., & Hughey J.J., 2021, PLoS Comput Biol
arguments
    Activity;
    binWidth = 0.5; % hour
    PeriodRange = [16 32]; % hour
    options.timePoints = (1:size(Activity,2)) * binWidth;
    options.ConfidenceLevel = 0.997; % 2 sigma ~ 95%; 3 sigma ~ 99.7%; 4 sigma ~ 99.994%; 5 sigma ~ 99.99994% 
    options.Method = 'Lomb-Scargle'; % 'Enright' ('Chi-Square'), 'Lomb-Scargle', 'Cosinor', 'greedy Chi-Square' 
end
options.Method = string(options.Method);
confidence = options.ConfidenceLevel;

P_candi = PeriodRange(1) : binWidth : PeriodRange(2); % candidate Periods for the periodogram
[rNum, N] = size(Activity); % row -> fly; column -> time points
t = options.timePoints;
powerArray = zeros(rNum, length(P_candi))*NaN;
peakPower = zeros(rNum,1)*NaN;
significanceFlag = zeros(rNum,1)*NaN;
pValue = zeros(rNum,1)*NaN;
period = zeros(rNum,1)*NaN;

if options.Method == "Enright" || options.Method == "Chi-Square"
    confidenceLine = chi2inv(confidence,P_candi/binWidth-1);
    for nn = 1:rNum
        M = mean(Activity(nn,:),'omitnan');
        for f = 1:length(P_candi)
            P_unit = P_candi(f)/binWidth; % period by the unit peirod (binWidth), 
            K = fix(N/P_unit);
            N_ef = K*P_unit;
            Mh = zeros(1,P_unit);
            for h = 1:P_unit
                Mh(h) = sum(Activity(nn,h:P_unit:(N_ef)),'omitnan') / K;
            end
            Qp = K * N_ef * ( sum((Mh-M).^2,'omitnan') ./ sum((Activity(nn,1:N_ef)-M).^2,'omitnan') );
            powerArray(nn,f) = Qp;
        end
        [m,I] = max(powerArray(nn,:) - confidenceLine);
        peakPower(nn) = powerArray(nn,I);
        pValue(nn) = 1 - chi2cdf(peakPower(nn),P_candi(I)/binWidth-1);
        if m > 0
            significanceFlag(nn) = 1;
            period(nn) = P_candi(I);
        else
            significanceFlag(nn) = 0;
        end
    end
    stat = struct();
    stat.Period = period;
    stat.PeakPower = peakPower;
    stat.SignificanceFlag = significanceFlag;
    stat.p_value= pValue;
    stat = orderfields(stat,{'Period','PeakPower','SignificanceFlag','p_value'});
elseif options.Method == "Lomb-Scargle"
    confidenceLine = repmat(-log(1-confidence^(1/N)), 1, length(P_candi));
    for nn = 1:rNum
        M = mean(Activity(nn,:),'omitnan');
        x = Activity(nn,:);
        for f = 1:length(P_candi)
            P = P_candi(f);
            delta = P/(4*pi) * atan(sum(sin(4*pi/P * t)) / sum(cos(4*pi/P * t)));
            A1 = (sum((x-M) .* cos(2*pi/P * (t-delta)))) ^ 2;
            B1 = sum((cos(2*pi/P * (t-delta))).^2);
            A2 = (sum((x-M) .* sin(2*pi/P * (t-delta)))) ^ 2;
            B2 = sum((sin(2*pi/P * (t-delta))).^2);
            PN = 1/(2*var(x)) * (A1/B1 + A2/B2);
            powerArray(nn,f) = PN;
        end
        [m,I] = max(powerArray(nn,:) - confidenceLine);
        peakPower(nn) = powerArray(nn,I);
        pValue(nn) = 1 - (1 - exp(-peakPower(nn)))^N;
        if m > 0
            significanceFlag(nn) = 1;
            period(nn) = P_candi(I);
        else
            significanceFlag(nn) = 0;
        end
    end
    stat = struct();
    stat.Period = period;
    stat.PeakPower = peakPower;
    stat.SignificanceFlag = significanceFlag;
    stat.p_value= pValue;
    stat = orderfields(stat,{'Period','PeakPower','SignificanceFlag','p_value'});
elseif options.Method == "Cosinor"
    confidenceLine = repmat(finv(confidence,2,N-3), 1, length(P_candi));
    A_array = zeros(rNum, length(P_candi))*NaN;
    phi_array = zeros(rNum, length(P_candi))*NaN;
    A_estimate = zeros(rNum,1)*NaN;
    phi_estimate = zeros(rNum,1)*NaN;
    for nn = 1:rNum
        Y_bar = mean(Activity(nn,:),'omitnan');
        Y = Activity(nn,:);
        for f = 1:length(P_candi)
            P = P_candi(f);
            x = cos(2*pi*t/P);
            z = sin(2*pi*t/P);
            d = [sum(Y);
                 sum(Y.*x);
                 sum(Y.*z)];
            S = [N,      sum(x),    sum(z);
                 sum(x), sum(x.^2), sum(x.*z);
                 sum(z), sum(x.*z), sum(z.^2)];
            u_est = S^(-1) * d; % u estimate
            
            M_est     = u_est(1);
            beta_est  = u_est(2);
            gamma_est = u_est(3);

            A_est = (beta_est^2 + gamma_est^2) ^ (1/2);
            phi_est = mod(atan(-gamma_est/beta_est),pi);
            Y_est = M_est + beta_est*x + gamma_est*z;
            RSS = sum((Y-Y_est).^2);
            MSS = sum((Y_est-Y_bar).^2);
            F = (MSS/2) / (RSS/(N-3));
            powerArray(nn,f) = F;
            A_array(nn,f) = A_est;
            phi_array(nn,f) = phi_est;
        end
        [m,I] = max(powerArray(nn,:) - confidenceLine);
        peakPower(nn) = powerArray(nn,I);
        A_estimate(nn) = A_array(nn,I);
        phi_estimate(nn) = phi_array(nn,I);
        pValue(nn) = 1 - fcdf(peakPower(nn),2,N-3);
        if m > 0
            significanceFlag(nn) = 1;
            period(nn) = P_candi(I);
        else
            significanceFlag(nn) = 0;
        end
    end
    stat = struct();
    stat.Period = period;
    stat.PeakPower = peakPower;
    stat.SignificanceFlag = significanceFlag;
    stat.p_value= pValue;
    stat.A_estimate = A_estimate;
    stat.phi_estimate = phi_estimate;
    stat.phi_by_pi = phi_estimate/pi;
    stat = orderfields(stat,{'Period','PeakPower','SignificanceFlag','p_value','A_estimate','phi_estimate','phi_by_pi'});
elseif options.Method == "greedy Chi-Square"
    confidenceLine = chi2inv(confidence,P_candi/binWidth-1);
    for nn = 1:rNum
        M = mean(Activity(nn,:),'omitnan');
        for f = 1:length(P_candi)
            P_unit = P_candi(f)/binWidth; % period by the unit peirod (binWidth), 
            K = N/P_unit;
            Mh = zeros(1,P_unit);
            for h = 1:P_unit
                Mh(h) = sum(Activity(nn,h:P_unit:(N)),'omitnan') / K;
            end
            Qp = K * N * ( sum((Mh-M).^2,'omitnan') ./ sum((Activity(nn,:)-M).^2,'omitnan') );
            powerArray(nn,f) = Qp;
        end
        [m,I] = max(powerArray(nn,:) - confidenceLine);
        peakPower(nn) = powerArray(nn,I);
        pValue(nn) = 1 - chi2cdf(peakPower(nn),P_candi(I)/binWidth-1);
        if m > 0
            significanceFlag(nn) = 1;
            period(nn) = P_candi(I);
        else
            significanceFlag(nn) = 0;
        end
    end
    stat = struct();
    stat.Period = period;
    stat.PeakPower = peakPower;
    stat.SignificanceFlag = significanceFlag;
    stat.p_value= pValue;
    stat = orderfields(stat,{'Period','PeakPower','SignificanceFlag','p_value'});
end

varargout = {stat, powerArray, confidenceLine, confidence, P_candi};
end




