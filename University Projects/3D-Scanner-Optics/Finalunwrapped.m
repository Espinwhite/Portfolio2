% 3-D Scanner by applying robust two-dimensional phase unwrapping
% Author: Alejandro Antonio Espinosa Gil
clc;  close all;

% Leemos las imagenes sin con objeto 
A1s = imread('0s.jpg'); A2s = imread('pi2s.jpg');
A3s = imread("pis.jpg"); A4s = imread("3pi2s.jpg");

% Leemos las imagenes con objeto 
A1 = imread('0.jpg'); A2 = imread('pi2.jpg');
A3 = imread("pi.jpg"); A4 = imread("3pi2.jpg");

% Convertir a blanco y negro con objeto
A1g = rgb2gray(A1); A2g = rgb2gray(A2);
A3g = rgb2gray(A3); A4g = rgb2gray(A4);

% Convertir a blanco y negro sin objeto
A1s = rgb2gray(A1s); A2s = rgb2gray(A2s);
A3s = rgb2gray(A3s); A4s = rgb2gray(A4s);

% Convertir a valores numéricos
I1s = double(A1s); I2s = double(A2s);
I3s = double(A3s); I4s = double(A4s);

% Convertir a valores numéricos
I1c = double(A1g); I2c = double(A2g);
I3c = double(A3g); I4c = double(A4g);

phi1 = atan2(I4c-I2c, I1c-I3c);    % Fase envuelta con objeto
phi2 = atan2(I4s-I2s, I1s-I3s);    % Fase envuelta sin objeto

%% Fase desenvuelta 
[y1] = fglobal(phi1); % Fase desenvuelta con objeto
[y2] = fglobal(phi2); % Fase desenvuelta sin objeto
phic1 = abs(y2-y1);                % Se resta la fase para obtener el objeto original 
phic2 = phic1';                    % Se genera otro phi para limpieza
m = size(phic2);

%% Limpieza Huevo
for i=1:m(1)
    for j=1:m(2)
        if (phic2(i,j)>3.3)
            phic2(i,j) =3.3;  
        elseif j < 170 || j > 1163 || i<314 || i>1756
            phic2(i,j) =0;
        else;continue;
        end
    end 
end

%limpieza tapa
for i=1:m(1)
    for j=1:m(2)
        if (phic2(i,j)>2)
            phic2(i,j) =0; 
        elseif j<=510 && phic2(i,j) >1
            phic2(i,j) = 0;
        elseif j>=640 && phic2(i,j) >1
            phic2(i,j) = 0;
        elseif i<=410 && phic2(i,j) >1
            phic2(i,j) = 0;
        elseif i>=770 && phic2(i,j) >1
            phic2(i,j) = 0;
        elseif phic2(i,j) <=0.5
            phic2(i,j) = 0;
        else;continue;
        end
    end 
end

%% Limpieza artefacto
for i=1:m(1)
    for j=1:m(2)
        if j>1015
            phic2(i,j)=0;

        elseif i>1142
            phic2(i,j)=0;
        elseif i<120
            phic2(i,j)=0;
        elseif (phic2(i,j)<=0.4)
                phic2(i,j)=0;
        else;continue;
        end
    end 
end

%% Limpieza muñeco
for i=1:m(1)
    for j=1:m(2)
        if (phic2(i,j)<=0.8)
           phic2(i,j)=0;
        elseif (phic2(i,j)>=3.9)
             phic2(i,j)=0;
        elseif j<760 && i>990
            phic2(i,j)=0;
        elseif j<760 && i<621
            phic2(i,j)=0;
        elseif j>900 && (phic2(i,j)>3.2)
            phic2(i,j)=0;
        elseif i>900 && (phic2(i,j)>3.2)
            phic2(i,j)=0;
        elseif i<260 
            phic2(i,j)=0;
        elseif i>1333
            phic2(i,j)=0;
        else;continue;
        end
    end 
end

%% Reescalamiento huevo
maxv = zeros(1,m(1));
for i=1:m(1)
    %phic2(i,:) = smooth(phic2(i,:),3); %Suamizamiento
    maxv(i) = max(phic2(i,:)); %Se calcula el maximo de cada vector  
end
maxg = max(maxv);
phic2 = phic2.*(4/maxg);
phic1 = phic1.*(4/maxg);

%% Reescalamiento muñeco
maxv = zeros(1,m(1));
for i=1:m(1)
    maxv(i) = max(phic2(i,:)); %Se calcula el maximo de cada vector  
end
maxg = max(maxv);
phic2 = phic2.*(1.9/maxg);
phic1 = phic1.*(1.9/maxg);

%% Smooth tapa
for i=1:m(1)
    phic2(i,:) = smooth(phic2(i,:),15); %Suamizamiento
end

%% Medidas originales huevo
xvec = linspace(0,15.5,m(2));
yvec = linspace(0,22.6,m(1));
zvec = linspace(-pi,pi,5);
[xvec1,yvec1] = meshgrid(xvec,yvec);

%% Medidas originales tapa
xvec = linspace(-6,13.6,m(2));
yvec = linspace(-6,14,m(1));
zvec = linspace(-pi,pi,5);
[xvec1,yvec1] = meshgrid(xvec,yvec);

%% Medidas originales plastilina muñeco
xvec = linspace(0,18.4,m(2));
yvec = linspace(0,10.9,m(1));
zvec = linspace(-pi,pi,5);
[xvec1,yvec1] = meshgrid(xvec,yvec);

%% Medidas originales plastilina artefacto
xvec = linspace(0,10.5,m(2));
yvec = linspace(0,8.5,m(1));
zvec = linspace(-pi,pi,5);
[xvec1,yvec1] = meshgrid(xvec,yvec);

%% Ejemplo de unwrapped huevo lineal
% Wrapped phase
figure(1)
plot(xvec,phi1(:,1020)), title('Fase envuelta'), xlabel('x (cm)'), ylabel('Radianes'), grid minor, yticks([-pi -pi/2 0 pi/2 pi]), yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% Unwrapped phase
figure(2)
plot(xvec,y1(:,1020)), title('Fase desenvuelta'), xlabel('x (cm)'), ylabel('Radianes'), grid minor,  yticks([0 pi/4 pi/2 3/4*pi pi 5/4*pi 3/2*pi]), yticklabels({'0','\pi/4','\pi/2','3/4\pi', '\pi','5/4\pi','3/2'\pi})
% Objeto real
figure(3)
plot(xvec,phic1(:,1020)), title('Fase desenvuelta'), xlabel('x (cm)'), ylabel('z (cm)'), grid minor

%% Ejemplo de unwrapped tapa lineal
figure(1)
plot(xvec,phi1(:,550)), title('Fase envuelta (tapa)'), xlabel('x (cm)'), ylabel('Radianes'),grid minor,  yticks([-pi -3/4*pi -pi/2 -pi/4 0 pi/2 pi]), yticklabels({'-\pi','-3/4\pi','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3/4\pi','\pi'})
figure(2)
plot(xvec,phic1(:,550)), title('Fase desenvuelta (tapa)'), xlabel('x (cm)'), ylabel('Radianes'), grid minor, yticks([0 pi/4 pi/2 3/4*pi pi 5/4*pi 3/2*pi]), yticklabels({'0','\pi/4','\pi/2','3/4\pi', '\pi','5/4\pi','3/2'\pi})

%% Ejemplo de unwrapped muñeco lineal
figure(1)
plot(xvec,phi1(:,800)), title('Fase envuelta (tapa)'), xlabel('x (cm)'), ylabel('Radianes'),grid minor,  yticks([-pi -3/4*pi -pi/2 -pi/4 0 pi/2 pi]), yticklabels({'-\pi','-3/4\pi','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3/4\pi','\pi'})
figure(2)
plot(xvec,phic1(:,800)), title('Fase desenvuelta (tapa)'), xlabel('x (cm)'), ylabel('Radianes'), grid minor, yticks([0 pi/4 pi/2 3/4*pi pi 5/4*pi 3/2*pi]), yticklabels({'0','\pi/4','\pi/2','3/4\pi', '\pi','5/4\pi','3/2'\pi})

%% Ejemplo de unwrapped artefacto
figure(1)
plot(xvec,phi1(:,800)), title('Fase envuelta'), xlabel('x (cm)'), ylabel('Radianes'),grid minor,  yticks([-pi -3/4*pi -pi/2 -pi/4 0 pi/2 pi]), yticklabels({'-\pi','-3/4\pi','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3/4\pi','\pi'})
figure(2)
plot(xvec,phic1(:,800)), title('Fase desenvuelta'), xlabel('x (cm)'), ylabel('Radianes'), grid minor, yticks([0 pi/4 pi/2 3/4*pi pi 5/4*pi 3/2*pi]), yticklabels({'0','\pi/4','\pi/2','3/4\pi', '\pi','5/4\pi','3/2\pi'})
figure(2)
plot(xvec,phic2(850,:)), title('Fase desenvuelta (tapa)'), xlabel('x (cm)'), ylabel('z [cm]'), grid minor

%% Surf de fase desnevuelta limpia tapa muñeco
figure(6)
surf(xvec1,yvec1,phic2), axis equal ,colormap white, shading interp; title('Fase desenvuelta Limpia'),xlabel('x [cm]'), ylabel('y [cm]'),zlabel('z[cm]'), clim([0 5]), camproj('perspective')
light("Style","local","Position",[0 0 4])

%% Surf de fase desnevuelta limpia tapa huevo
figure(6)
surf(xvec1,yvec1,phic2), axis equal ,colormap white, shading interp; title('Fase desenvuelta Limpia'),xlabel('x [cm]'), ylabel('y [cm]'),zlabel('z[cm]'), clim([0 5]), camproj('perspective')
light("Style","local","Position",[0 0 4])

%% Surf fase desenvbuelta tapa
y = linspace(0,9,9); x = linspace(0,8,8)
figure(6)
surf(xvec1,yvec1,phic2), axis equal ,colormap white, shading interp; title('Fase desenvuelta Limpia'),xlabel('x [cm]'), ylabel('y [cm]'),zlabel('z(cm)'), clim([0 5]), camproj('perspective'),xlim([0 8]), ylim([0 9])
light("Style","local","Position",[0 0 4])
toc

%% Surf de fase desnevuelta limpia artefacto
figure(6)
surf(xvec1,yvec1,phic2), axis equal ,colormap white, shading interp; title('Fase desenvuelta Limpia'),xlabel('x'), ylabel('y'),zlabel('z(cm)'), clim([0 5]), camproj('perspective'), xlim([0 8]), ylim([0 8])
light("Style","local","Position",[0 0 4])

function [phaseu] = fglobal(phasew) 
   phii = finterna(phasew);           %Mandamos a llamar la función para estimar el primer desenvolvimiento de la fase (phi1)
   phii = phii + mean2(phasew) - mean2(phii); %mean2 calcula la media de todos los valores de la primera etapa desenvuelta 
   k1 = round((phii-phasew)/2/pi);        %Se calcula el entero k1
   phaseu = phasew + 2*k1*pi;         %Se utiliza el entero k1 para calcular la fase desenvuelta
   error = wrapToPi(phaseu-phii);     %WraptoPi te acota tus fases entre un rango de [-pi pi]
   phii = phii + finterna(error);           %Se evalua el error
   phii = phii + mean2(phasew)-mean2(phii); %Se actualiza phii
   k2 = round((phii-phasew)/2/pi);        %Se calcula el entero k2
   phaseu = phasew + 2*k2*pi;         %Se contruye de nuevo la fase desenvuelta con el entero k2
   error = wrapToPi(phaseu-phii);     %WraptoPi te acota tus fases entre un rango de [-pi pi]
   n = 0;                                     %Se actualiza el numero de iteraciones
   while sum(sum(abs(k2-k1)))>0             %Condición para dejar de iterar
    k1 = k2;                                  %Se actualiza k1
    phic = finterna(error);               %Se aplica de nuevo la funcnión interna
    phii = phii + phic;                         %Se actualiza la fase, utilizanod la fase correcta
    phii = phii + mean2(phasew)-mean2(phii);%adjust piston
    k2 = round((phii-phasew)/2/pi);       %Se calcula la segunda K
    phaseu = phasew + 2*k2*pi;        %Se calcula la segunda fase desenvuelta de nuestro algoritmo recursivo
    error = wrapToPi(phaseu-phii);    %Se calcula el error
    n = n + 1;                                  %Se actualiza el numero de iteraciones
   end
end
% Función para desenvolver ecuación de transporte de intensidad 
function [phaseu] = finterna(phasew)
    psiw = exp(1i*phasew);               %Se contruye la exponente compleja con la fase envuelta
    edx = [zeros([size(psiw,1),1]), wrapToPi(diff(psiw, 1, 2)), zeros([size(psiw,1),1])]; %Se genera una función de x y y para la fase envuelta de la exponente compleja
    edy = [zeros([1,size(psiw,2)]); wrapToPi(diff(psiw, 1, 1)); zeros([1,size(psiw,2)])];
    den = diff(edx, 1, 2) + diff(edy, 1, 1); %Se calcula el laplaciano de la exponente con la fase envuelta
    rho = imag(conj(psiw).*den);             % Se calcula el numerador de la ecuación (9)
    phaseu = junta(rho);       %Se computa la ecuación (9) mandando llamar a solvePoisson, para juntar el numerador y denominador 
end
%La funcnión "junta" junta el numerador y el denominador de la
%ecuación (9) para obtener la fase desenvuleta
function phii = junta(rho)
    cphi2 = dct2(rho);                     % Resuleve la ecuación de poisson con tranformación de coseno discreta % en 2d
    [N, M] = size(rho);                     %Se generan los vectores "N" y "M"
    [I, J] = meshgrid(0:M-1, 0:N-1);    %vectores "x" y "y" 
    cPhi2 = cphi2 ./ 2 ./ (cos(pi*I/M) + cos(pi*J/N) - 2); %Se computa la ecuación (9) del paper
    cPhi2(1,1) = 0;                        % Se renombra la posición (1,1) de la trasnfromación de coseno discreta
    phii = idct2(cPhi2);                    % Se obtiene la tranformación de coseno discreta inversa 
end