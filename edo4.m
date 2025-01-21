function dudt=edo4m(tinc2,u,m,k,C,incl1)
dudt = [u(2); (C*(-incl1)+k*(-incl1*(tinc2-1.8)) - C*u(2) - k*u(1))/m];
%a forca e a forca total de reacao da pista que esta agindo sobre o
%carrinho
