function [Loc] = mat_local(Nw, Nm, obs, L, Lx, Ly, theta)

n = sqrt(Nm);
Loc = zeros(Nm, Nw);

trans = zeros(2);
trans(1, 1) = cos(theta); trans(1, 2) = -sin(theta);
trans(2, 1) = sin(theta); trans(2, 2) = cos(theta);

for kk = 1:Nw
    
    Loc_aux = zeros(Nm, 1);
    cont = 1;
    
    for i = 1:n
        for j = 1:n
            
            auxi = (i - obs(kk, 1));
            auxj = (j - obs(kk, 2));
            vec = trans*[auxi; auxj];
            h = sqrt((vec(1)/Lx)^2 + (vec(2)/Ly)^2);
            Loc_aux(cont) = rho(h, L);
            cont = cont + 1;
        end
    end
    Loc(:, kk) = Loc_aux;
end


end

function [r] = rho(h, L)

if (and( h >= 0, h <= L))
    r = -(1/4)*(h/L)^5 + (1/2)*(h/L)^4 + (5/8)*(h/L)^3 - (5/3)*(h/L)^2 + 1;
elseif (and(h > L, h <= 2*L))
    r = (1/12)*(h/L)^5 - (1/2)*(h/L)^4 + (5/8)*(h/L)^3 + (5/3)*(h/L)^2 - 5*(h/L) + 4 - (2/3)*(h/L)^(-1);
else
    r = 0;
end

end