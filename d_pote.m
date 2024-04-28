function dP = d_pote(d, dD, dM, lam, mu)
    if d <= dD
%         dP = (d - dD)*(d + dD + 2*dD^2/mu)/(d + dD^2/mu)^2;
        dP = -(dD - d)^2*(2*d + dD + 3*dD^3/mu)/(d + dD^3/mu)^2;
    elseif  dD < d && d <= dM
        dP = lam*d_smstep(d, dD, dM);
    else
        dP = 0;
    end
end