clc; clear; close all;

n = 6; % dimensione e numero di matrici Qk
[x_opt,Qlist,Qprod] = RoundRobinQmatrix(n);

disp('Soluzione ottimale per le variabili [x1..x2n]:')
disp(x_opt.')

for k = 1:n
    fprintf('Matrice Q%d:\n',k)
    disp(Qlist{k})
end

disp(['Prodotto (Q1*...*Qn)^2 (Check convergenza)'])
disp(['NB. mi aspetto di ottenere una matrice con tutti i termini pari a: ',num2str(1/n)])
disp(Qprod^2)

disp('Somma righe di Qprod:')
disp(sum(Qprod,2).')

disp('Somma colonne di Qprod:')
disp(sum(Qprod,1))