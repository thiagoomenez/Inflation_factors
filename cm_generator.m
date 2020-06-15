function [c] = cm_generator(nx, ny, theta, Lmax, Lmin, deltak)

n = nx*ny;
c = zeros(n);

transfer = zeros(2);
transfer(1, 1) = cos(theta);
transfer(1, 2) = sin(theta);
transfer(2, 1) = - sin(theta);
transfer(2, 2) = cos(theta);

for i = 1:n
    for j = i:n
        [vaux1i, vaux1j] = cm_find(i, nx, ny);
        [vaux2i, vaux2j] = cm_find(j, nx, ny);
        vaux = [(vaux1i - vaux2i); (vaux1j - vaux2j)];
        v = transfer*vaux;
        h = sqrt( (v(1)/Lmax)^2 + (v(2)/Lmin)^2 );
        cov = spherecov(h);
        c(i, j) = cov*deltak*deltak;
        c(j, i) = cov*deltak*deltak;
    end
end

c = chol(c, 'lower');

end

function [a, b] = cm_find(x, nx, ny)

cont = 1;
caux = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        caux(i, j) = cont;
        cont = cont+1;
    end
end

u = find(caux == x);

a = mod(u, nx);
if (a == 0)
    a = nx;
end
b = ceil(u/ny);

end

function [s] = spherecov(h)

    if(h < 1)
        s = (1 - (3/2)*h + h*h*h/2);
    else
        s = 0;
    end

end