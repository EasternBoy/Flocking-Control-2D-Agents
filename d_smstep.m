function y = d_smstep(x,a,b)    
    z = (x-a)/(b-a);
    if x < a
        y = 0;
    elseif x > b
        y = 0;
    else
        y = 30*z^4 - 60*z^3 + 30*z^2;
    end
end