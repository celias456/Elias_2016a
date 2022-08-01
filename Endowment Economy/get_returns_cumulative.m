function [cum_1_percent,cum_2_percent,cum_4_percent] = get_returns_cumulative(gross)
%Calculate 1,2,and 4 year ahead cumulative return (percentage) from gross returns

%Accepts as input a series of gross returns

n = length(gross);

cum_1 = zeros(n,1);
cum_2 = zeros(n,1);
cum_4 = zeros(n,1);

for i = 1:n
    if i >= 4
        cum_1(i) = gross(i)*gross(i-1)*gross(i-2)*gross(i-3);
    end
    if i >= 8
        cum_2(i) = gross(i)*gross(i-1)*gross(i-2)*gross(i-3)*gross(i-4)*gross(i-5)*gross(i-6)*gross(i-7);
    end
    if i >= 16
        cum_4(i) = gross(i)*gross(i-1)*gross(i-2)*gross(i-3)*gross(i-4)*gross(i-5)*gross(i-6)*gross(i-7)*gross(i-8)*gross(i-9)*gross(i-10)*gross(i-11)*gross(i-12)*gross(i-13)*gross(i-14)*gross(i-15);
    end   
            
end

cum_1 = cum_1(4:end);
cum_2 = cum_2(8:end);
cum_4 = cum_4(16:end);

cum_1_percent = (cum_1 - 1)*100;
cum_2_percent = (cum_2 - 1)*100;
cum_4_percent = (cum_4 - 1)*100;


end

