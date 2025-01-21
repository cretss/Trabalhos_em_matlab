function dudt=EDO2m(u,m,k,C)
dudt = [u(2); (- C*(u(2)) - k*(u(1)))/m];
end
