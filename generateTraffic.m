function [tlambda] = generateTraffic(arrival, duration, cv, lambda)

if strcmp(arrival,'M')
    tlambda = exprnd(1/lambda,1,duration); % s
elseif strcmp(arrival,'G')
    m=1/lambda;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tlambda = lognrnd(mx,sigma,1,duration);
elseif strcmp(arrival,'AR')
    m=1/lambda;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tlambda = lognrnd(mx,sigma,1,duration);
    for i = 2:duration
        tlambda(i) = tlambda(i-1) * autocorrelationvalue + tlambda(i) * (1-autocorrelationvalue);
    end
end
