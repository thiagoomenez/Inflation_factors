function [m_a] = esmda(Ne, m_f, d, alpha, Cd, dobs, Loc)

Cmd = (1/(Ne-1))*(m_f - mean(m_f, 2))*(d - mean(d, 2))';
Cdd = (1/(Ne-1))*(d - mean(d, 2))*(d - mean(d, 2))';
K = Cmd/(Cdd + alpha*Cd);

m_a = m_f + (Loc.*K)*(dobs + sqrt(alpha*Cd)*randn(Nd, 1) - d);

end
