function [tmu] = generateVideo(service, duration, cv, mu)

if strcmp(service,'M')
    tmu = exprnd(mu,1,duration); % s
elseif strcmp(service,'G')
    m=mu;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tmu = lognrnd(mx,sigma,1,duration);
elseif strcmp(service,'AR')
    m=mu;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tmu = lognrnd(mx,sigma,1,duration);
    for i = 2:duration
        tmu(i) = tmu(i-1) * autocorrelationvalue + tmu(i) * (1-autocorrelationvalue);
    end
end
