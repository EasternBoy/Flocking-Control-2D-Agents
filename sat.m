function out = sat(in,para)
    temp = zeros(max(size(in)),1);
    for i = 1:1:max(size(in))
        temp(i) = in(i)/sqrt(para^2 + in(i)^2);
    end
    out = temp;
end

