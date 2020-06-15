function [m] = esmda(Ne, m, d, alpha, Cd, dobs, Loc)

Cmd = (1/(Ne-1))*(m - mean(m, 2))*(d - mean(d, 2))';
Cdd = (1/(Ne-1))*(d - mean(d, 2))*(d - mean(d, 2))';
K = Cmd/(Cdd + alpha*Cd);

m = m + (Loc.*K)*(dobs + sqrt(alpha*Cd)*randn(Nd, 1) - d);

end