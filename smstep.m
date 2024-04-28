function y = smstep(x,a,b)   

    z = (x-a)/(b-a);
    if x <= a
        y = 0;
    elseif x >= b
        y = 1;
    else
        y = 6*z^5 - 15*z^4 + 10*z^3;
    end
end