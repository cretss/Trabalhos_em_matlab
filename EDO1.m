function dudt=EDO1m(tn,u,m,k,C,f,A)
dudt = [u(2); (C*(A/2)*2*pi*f*cos(2*pi*f*tn)+k*(A/2)*sin(2*pi*f*tn) -
C*u(2) - k*u(1))/m];
%a forca e a forca total de reacao da pista que esta agindo sobre o
%carrinho
