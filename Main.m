%Parâmetros do sistema
t= 0:0.1:10; %s
lomb=25;%comprimento da lombada
A=0.08; %m (Amplitude de 8 cm)
S=50; %m (comprimento da onda)
v=11.11; % m/s Velocidade do carro (40Km/h)
f=(v/S); %Hz,
k=80000; %N/m
c=5000; %N.s/m
m=1000; %Kg
temp1=(1/f)/2;%tempo de percorrer a lombada 2
ta=0:0.01:temp1;
ts=temp1:0.01:10;
%achando frequencia natural do meu sistema e a de amortecimento
Wn=sqrt(k/m)
zeta=(c)/(2*sqrt(k*m))
Wd=2*pi*f;
r=Wn/Wd;
%frequencia harmonica da base com a ideia de variacao do r
rvar=0:0.01:10;
RT=sqrt(1+(2*zeta*rvar).^2)./sqrt((1-rvar.^2).^2+(2*zeta*rvar).^2);
figure(1)
plot(rvar,RT)
title('FRF do sistema por forçamento harmônico de base');
xlabel('Razão de frequências (wd/wn)');
ylabel('Razão de transmissibilidade');
grid on;
%Resposta da massa como um impulso e seu deslocamento
x0=0; %m
dx0=4; %m/s
x=exp(-zeta*Wn.*t).*(x0*cos(Wd.*t)+(dx0+zeta*Wn*x0)./(Wd)*sin(Wd.*t));
figure(2)
plot(t,x)
title('Resposta Impulsiva do sistema');
xlabel('Tempo(s)');
ylabel('Deslocamento vertical da massa u(m) sujeita a impulso');
grid on;
%analise da pista primeira pista
%Inicializando os vetores necessários:
%Equação para primeira pista:
y1=(A/2)*sin(2*pi*f*t);
dy1=(A/2)*2*pi*f*cos(2*pi*f*t);
%resolucao da edo
[t1,U]=ode45(@(tn,u) EDO1m(tn,u,m,k,c,f,A),0:0.1:10,[0;0]);
%Plotando
figure(3)
plot(t,y1)
ylim([-0.1 0.1])
title('Modelagem da pista 1');
xlabel('Tempo (s)');
ylabel('y(t)');
grid on;
figure(4)
plot(t,dy1);
title('Derivada da modelagem da pista 1');
xlabel('Tempo (s)');
ylabel('dy(t)');
grid on;
figure(5)
plot(t1,U(:,1));
ylim([-0.06 0.06])
title('Movimento da massa na pista1')
xlabel('Tempo (s)'); ylabel('Deslocamento (m)');
grid on;
%pista 2
%criando parametros pra isso
y2=(A/2)*sin(2*pi*f*ta);
dy2=(A/2)*2*pi*f*cos(2*pi*f*ta);
y3=zeros(1,numel(ts),'double');
dy3=zeros(1,numel(ts),'double');
figure(6)
plot(ta,y2,ts,y3)
title('Modelagem da pista 2');
ylim([-0.1 0.1])
xlabel('Tempo (s)');
ylabel('y(t)');
grid on;
figure(7)
plot(ta,dy2,ts,dy3);
ylim([-20 20])
title('Derivada da Modelagem da pista 2 ');
xlabel('Tempo (s)');
ylabel('dy(t)');
grid on;
%Resolução da EDO
[t12,U2]=ode45(@(tn,u) EDO1m(tn,u,m,k,c,f,A),0:0.01:temp1,[0;0]);
u0=9.14055e-5;%condicoes de entrada da edo2
u01=-0.0573783;%condicoes de entrada da edo2
[t13,U3]=ode45(@(tn,u) EDO2m(u,m,k,c),temp1:0.01:10,[u0;u01]);
figure(8)
plot(t12,U2(:,1),t13,U3(:,1));
ylim([-0.02 0.06])
title('Movimento da massa na pista 2')
xlabel('Tempo (s)'); ylabel('Deslocamento (m)');
grid on;
%figure(8)grafico para achar as condicoes de contorno
%plot(U2(:,2),U2(:,1));
%ylim([-0.02 0.06])
%title('Movimento da massa na pista 2')
%xlabel('Tempo (s)'); ylabel('Deslocamento (m)');
%grid on;
%%%%pensando na lombada trapezoidal
%%%precisa modelar o trapezio
%25 total -5 parte de subida e descida 1.5 tamanho
incl1=A/(5/11.11);
tinc=0:0.01:0.45;
tficar=0.45:0.01:1.8;
tinc2=1.8:0.01:2.25;
tend=2.25:0.01:5;
y4=incl1*tinc;
dy4=incl1*ones(1,numel(tinc),'double');
y5=(incl1*0.45)*ones(1,numel(tficar),'double');
dy5=zeros(1,numel(tficar),'double');
y6=(incl1*0.45)-incl1*(tinc2-1.8);
dy6=-incl1*ones(1,numel(tinc2),'double');
y7=zeros(1,numel(tend),'double');
dy7=zeros(1,numel(tend),'double');
figure(9)
plot(tinc,y4,tficar,y5,tinc2,y6,tend,y7)
title('Modelagem da pista 3');
ylim([-0.1 0.1])
xlabel('Tempo (s)');
ylabel('y(t)');
grid on;
figure(10)
plot(tinc,dy4,tficar,dy5,tinc2,dy6,tend,dy7)
ylim([-0.4 0.4])
title('Derivada da Modelagem da pista 3');
xlabel('Tempo (s)');
ylabel('dy(t)');
grid on;
[t14,U4]=ode45(@(tn,u) edo3m(tn,u,m,k,c,incl1),0:0.01:0.45,[0;0]);
u0s=0.0844377;
u01s=0.209922;
[t15,U5]=ode45(@(tn,u) EDO6m(u,m,k,c),0.45:0.01:1.8,[u0s;u01s]);
u021s=0.000692026;
u021=0.0793585;
[t16,U6]=ode45(@(tn,u) edo7m(tn,u,m,k,c,incl1),1.8:0.01:2.25,[u021;u021s]);
u03s=-0.214512;
u03=-0.0852102;
u031s=-0.211331;
u031=-0.0042587;
%[t17,U7]=ode45(@(tn,u) EDO2m(u,m,k,c),2.25:0.01:5,[u03;u03s]);
[t17,U7]=ode45(@(tn,u) EDO2m(u,m,k,c),2.25:0.01:5,[u031;u031s]);
figure(11)
%plot(t16,U6(:,1),t16,U6(:,2)) %para achar as condicoes iniciais
plot(t14,U4(:,1),t15,U5(:,1),t16,U6(:,1),t17,U7(:,1));
ylim([-0.15 0.15])
title('Movimento da massa na pista 3')
xlabel('Tempo (s)'); ylabel('Deslocamento (m)');
grid on;
