function P = pote(d, dD, dM, lam, mu)
%    P = (d-dD)^2*(1 - smstep(d, dD, dM))/d^2 + lam*smstep(d, dD, dM);
%    P =  (d-dD)^2*(1 - smstep(d, dD, dM))/d + lam*smstep(d, dD, dM);
    if d <= dD
        P = (d-dD)^2/(d + dD^2/mu);
%         P = -(d-dD)^3/(d + dD^3/mu);
%         P = (d-dD)^2/d^2;
    elseif d >= dM
        P = lam;
    else
        P = smstep(d, dD, dM);
    end
end