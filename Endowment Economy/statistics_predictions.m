function [error] = statistics_predictions(allvariables_levels)
%Calculates prediction statistics

price = allvariables_levels(:,10); %previous period price
predict1 = allvariables_levels(:,11); %Agent 1 previous period prediction
predict2 = allvariables_levels(:,12); %Agent 2 previous period prediction
n = length(price);

error1 = (price - predict1).^2;
error1_sum = sum(error1);
mse_1 = error1_sum/n;

error2 = (price - predict2).^2;
error2_sum = sum(error2);
mse_2 = error2_sum/n;

error = [mse_1, mse_2];



end

