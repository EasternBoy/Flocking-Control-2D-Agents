function [dij, Gij, Hij] = GH_ij(qi,qj,qod,poi,poj,Mi,Mj,dD,dM,lam,mu)

    phi_i = qi(3);
    phi_j = qj(3);
    pi    = qi(1:2);
    pj    = qj(1:2);
    pod   = qod(1:2);
    Ri    = [cos(phi_i) -sin(phi_i);sin(phi_i) cos(phi_i)];
    Rj    = [cos(phi_j) -sin(phi_j);sin(phi_j) cos(phi_j)];
    S     = [0 -1; 1 0];

    Den_ij = 0;
    
    delij = zeros(Mi,Mj,Mj);
    delji = zeros(Mj,Mi,Mi);

%% Calculate relative distance dij
    for k = 1:Mi
        pki = pi + Ri*poi(:,k);
        for n = 1:Mj
            pnj = pj + Rj*poj(:,n);
            if n < Mj
                m = n + 1;
                pmj = pj + Rj*poj(:,n+1);
            else
                m = 1;
                pmj = pj + Rj*poj(:,1);
            end
            delij(k,n,m) = normE(pki-pnj)+normE(pki-pmj)-normE(pnj-pmj);
            Den_ij = Den_ij + 1/delij(k,n,m);       
        end
    end

    for k = 1:Mj
        pkj = pj +Rj*poj(:,k);
        for n = 1:Mi
            pni = pi + Ri*poi(:,n);
            if n < Mi
                m = n + 1;
                pmi = pi + Ri*poi(:,n+1);
            else
                m = 1;
                pmi = pi + Ri*poi(:,1);
            end
            delji(k,n,m) = normE(pkj-pni)+normE(pkj-pmi)-normE(pni-pmi);
            Den_ij = Den_ij + 1/delji(k,n,m);
        end
    end

    dij = 2*Mi*Mj*1/Den_ij;

%% Calculate Gij and Hij
if dij < dM  % Out of communication range

    Eij = [0 0];
    Fij = 0;
    cof = dij^2/(2*Mi*Mj);

    for k = 1:Mi
        pki = pi + Ri*poi(:,k);
        for n = 1:Mj
            pnj = pj + Rj*poj(:,n);
            if n < Mj
                m = n + 1;
                pmj = pj + Rj*poj(:,n+1);
            else
                m = 1;
                pmj = pj + Rj*poj(:,1);
            end
            Akmn = (pki-pnj)'/normE(pki-pnj) + (pki-pmj)'/normE(pki-pmj);
            temp = Akmn/delij(k,n,m)^2;
            Eij  = Eij + temp;
            Fij  = Fij + temp*S*(pki - pi);
        end
    end

    Eij = cof*Eij;
    Fij = cof*Fij;

    Eji = [0 0];
    Fji = 0;

    for k = 1:Mj
        pkj = pj + Rj*poj(:,k);
        for n = 1:Mi
            pni = pi + Ri*poi(:,n);
            if n < Mi
                m = n + 1;
                pmi = pi + Ri*poi(:,n+1);
            else
                m = 1;
                pmi = pi + Ri*poi(:,1);
            end
            Akmn = (pkj-pni)'/normE(pkj-pni) + (pkj-pmi)'/normE(pkj-pmi);
            temp = Akmn/delji(k,n,m)^2;
            Eji  = Eji + temp;
            Fji  = Fji + temp*S*(pkj - pj);
        end
    end
    
    Eji = cof*Eji;
    Fji = cof*Fji;

    Gij = d_pote(dij,dD,dM,lam,mu)*(Eij - Eji);
    Hij = d_pote(dij,dD,dM,lam,mu)*(Eij*S*(pi - pod) + Fij - Eji*S*(pj - pod) - Fji);
else
    Gij = [0 0];
    Hij = 0;
end

end

