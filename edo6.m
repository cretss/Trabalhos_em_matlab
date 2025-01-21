function dudt=EDO6m(u,m,k,C)
dudt = [u(2); (- C*(u(2)) - k*(u(1)-0.08))/m];
end
