clear all

global d1
global d2
global cB1
global cB2

charge = 1.66*10^-19;
kB = 1.38*10^-23;
T = 300;
eps = 78*8.85*10^-12;
Na = 6.022*10^23;
c_bar = 10^3 * Na;

lD = (eps*kB*T/2/charge^2/c_bar)^0.5;
DR=6.5;
R_rod = (DR/2)*10^-9/lD;
D = 9.6*10^-9/lD;
D1=D*sqrt(3)/2;


cB1 = c_bar * lD^3;
cB2 = c_bar * lD^3;
d1 = 3.82*2*10^-10/lD ;
d2 = 3.32*2*10^-10/lD;


V_rod = 0.2 * charge / kB / T;

global V_list
global Sol1List
global Sol2List

V_list = linspace(-3*V_rod,3*V_rod,2000); 
Sol1List = zeros(1, length(V_list));
Sol2List = zeros(1, length(V_list));

for i = 1:length(V_list)
    fun = @(t) root2d(t, V_list(i));
    t0 = [cB1, cB2];
    options = optimoptions('fsolve', 'Display', 'off');
    [X,fval,exitflag] = fsolve(fun, t0,options);
    Sol1List(i) = X(1);
    Sol2List(i) = X(2);
end

C2 = [1;0;0;5*R_rod]; 
C1 = [1;0;0;R_rod];
C3 = [1;0;D;R_rod];
C4 = [1;0;-D;R_rod];
C5=[1;D1;-D/2;R_rod];
C6=[1;D1;D/2;R_rod];
C7=[1;-D1;D/2;R_rod];
C8=[1;-D1;-D/2;R_rod];
gd = [C2,C1,C3,C4,C5,C6,C7,C8];
sf = 'C2-C1-C3-C4-C5-C6-C7-C8';
ns = char('C2','C1','C3','C4','C5','C6','C7','C8')';
g = decsg(gd,sf,ns);
model=createpde %создаем пустую модель
geometryFromEdges(model,g); %объединяем пустую модель и нашу геометрию
applyBoundaryCondition(model,'dirichlet','edge',[1:4],'u',0);
applyBoundaryCondition(model,'dirichlet','edge',[5:33],'u',V_rod); 
specifyCoefficients(model,'m',0,'d',0,'c',-1,'a',0,'f',@fcoeffunction); 
mesh=generateMesh(model,'Hmax',0.5) %генерируем Мэш
results = solvepde(model); %Решаем наше уравнение
u = results.NodalSolution; %Получаем численный результат
xq = linspace(R_rod,1.5*R_rod,100); %задаем х на которых будем интерполировать наше решение
nt= length(xq);
x=xq';
yq = zeros(1,nt); %интерполируем только вдоль оси х
uintr=interpolateSolution(results,xq,yq); %интерполяция
xy=horzcat(x,uintr);
sum=0;
listx=[];
i=1;
listy=[];
j=1;
deltay=0.1;
deltax=0.1;
for x1= -5*R_rod:deltax:5*R_rod
    for y1= -5*R_rod:deltay:5*R_rod
        k=interpolateSolution(results,x1,y1);
        if (~isnan(k))
             X_1 = interp1(V_list,Sol1List,k);
             X_2 = interp1(V_list,Sol2List,k);
             q=(X_1*charge-charge*X_2);
             sum=sum+q;
   
        end
        
  
    end
end
sum=sum*deltay*deltax;
S=7*2*pi*R_rod;
save('LiCl_D65_02_plusepsilon2.mat')
function f = fcoeffunction(location,state) %функция для коэффициента f
global cB1

global V_list
global Sol1List
global Sol2List

N = 1; 
nr = length(location.x);
f = zeros(N,nr);

for i = 1:nr
    if ( ~isnan( state.u(1,i) ) )
        X_1 = interp1(V_list,Sol1List,state.u(1,i));
        X_2 = interp1(V_list,Sol2List,state.u(1,i));
        f(1,i) = 0.5 * (X_2 - X_1)/cB1;
    else
        f(1,i) = NaN;
    end
end

end


function out = pR(c1,c2)
global d1
global d2
n0 = c1 + c2;
n1 = 0.5*d1*c1 + 0.5*d2*c2;
n2 = pi*d1^2*c1 + pi*d2^2*c2;
n3 = pi/6*d1^3*c1 + pi/6*d2^3*c2;
first = n0*(1 + n3 + n3^2)/(1-n3)^3;
second = (1/12/pi*n2^3 + n1*n2*(1-n3) - 3*n0*n3)/(1-n3)^3;
out = (first + second);
end

function out = chpRex_1(c1,c2)
global d1
global d2
n1 = 0.5*d1*c1 + 0.5*d2*c2;
n2 = pi*d1^2*c1 + pi*d2^2*c2;
n3 = pi/6*d1^3*c1 + pi/6*d2^3*c2;
first = -log(1 - n3);
second = n2*d1/2/(1 - n3);
third = ( n1/(1 - n3) + n2^2/8/pi/(1 - n3)^2 )*pi*d1^2;
fourth = pR(c1,c2)*pi/6*d1^3;
out = (first + second + third + fourth);
end

function out = chpRex_2(c1,c2)
global d1
global d2
n1 = 0.5*d1*c1 + 0.5*d2*c2;
n2 = pi*d1^2*c1 + pi*d2^2*c2;
n3 = pi/6*d1^3*c1 + pi/6*d2^3*c2;
first = -log(1 - n3);
second = n2*d2/2/(1 - n3);
third = ( n1/(1 - n3) + n2^2/8/pi/(1 - n3)^2 )*pi*d2^2;
fourth = pR(c1,c2)*pi/6*d2^3;
out = (first + second + third + fourth);
end

function F = root2d(x, psi)
global cB1
global cB2
F(1) = x(1) - cB1*exp(-chpRex_1(x(1),x(2)) + chpRex_1(cB1,cB2) - psi);
F(2) = x(2) - cB2*exp(-chpRex_2(x(1),x(2)) + chpRex_2(cB1,cB2) + psi);
end