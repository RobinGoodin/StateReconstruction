tic
clear
clc
t=100000; 
p=1/sqrt(2);
NUM = 1000;

Tau = 0:0.1:2;
T2 = 1;
result = zeros(length(Tau),1);

psi=[p;1i*p];
psi=psi/norm(psi);

po = psi*psi';
% Состояние задано
X=[1,0;0,1;p,p;p,-p;p,-1i*p;p,1i*p]; %X-Аппаратная матрица
[N N0]=size(X); % N-число базисов в аппаратной матрице

phi = zeros(N,1);
theta = zeros(N,1);
s = zeros(N,1);
for k=1:N
    phi(k) = phase(X(k,2))-phase(X(k,1));
	theta(k) = 2*acos(abs(X(k,1)));
    s(k) = (phase(X(k,1))+phase(X(k,2)))/2;
end

for QT = 1:length(Tau)
    
T = Tau(QT);
L_phase_T = [];

for k=1:N
    buffer = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k))*exp(-T/T2);
              sin(theta(k))*exp(-1i*phi(k))*exp(-T/T2),1-cos(theta(k))];
    buffer = buffer*1/2;
    L_phase_T(:,:,k) = buffer;
    %X(k,:)'*X(k,:) проверка для операторов измерения
end

L = zeros(N,1);
for k=1:N
    L(k,1)=trace(L_phase_T(:,:,k)*po);
end

res = zeros(NUM,1);
for QW = 1:NUM
% Генерация числа фотонов, подчиняющихся распределению Пуассона
K=zeros(N,1);
for k=1:N
    K(k,1)=poissrnd(L(k,1)*t);
end

% Неизменная матрица для нахождения С
I=zeros(N0,N0); 
I = I + X'*X*t;
% Решение методом итераций
% po1 - Первое приближение
po1 = [1 0; 0 1];
po1 = po1/trace(po1);
% a - Параметр сходимости
a=0.5; 
for k=1:1000
    L2=zeros(N,1);
    for k=1:N
        L2(k,1)=trace(L_phase_T(:,:,k)*po1);
    end
    % Матрица, зависящая от С
    J=zeros(N0,N0); 
    for j=1:N0
        for n=1:N0
            for m=1:N
                 J(j,n)=J(j,n)+conj(X(m,j))*X(m,n)*K(m)/L2(m);
            end
        end
    end
    po1=(1-a)*I^(-1)*J*po1+a*po1;
    po1=po1/trace(po1);
end
F = (trace((po^(1/2)*po1*po^(1/2))^(1/2)))^2;
res(QW) = real(1-F);

end

result(QT) = sum(res)/NUM;
QT
end
result
plot(Tau,result);

% po
% po1
% F=vpa((trace((po^(1/2)*po1*po^(1/2))^(1/2)))^2,10)

toc