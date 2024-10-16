close all
clear all



% DEFINITION DES PARAMETRES DU PROBLEME
params.fonction ='rosenbrock';
% A COMPLETER 
params.n= 2;
params.b = 2;
params.forme = 'moindres carrés';

% DEFINITION DE L ALGORITHME D OPTIMISATION 

options.pas ='variable';
% A COMPLETER
options.const = 10^-4;
options.beta = 0.75;
options.tolX = 10^-8;
options.tolF = 10^-8;
options.tolG = 10^-8;
options.maxiter = 10^3;

% OPTIMISATION

%AFFICHACHE
% Définition des méthodes et des couleurs

%methods = {'gradient', 'gradient conjuge', 'Newton', 'Quasi-Newton'};
methods = {'Gauss-Newton', 'Levenberg-Marquardt'};
colors = {'b', 'r', 'g', 'm'};  % Couleurs différentes pour chaque méthode

% Variables pour stocker les résultats
num_iterations = zeros(length(methods), 1);
execution_time = zeros(length(methods), 1);

figure;  % Créer une seule figure pour toutes les méthodes
hold on;  % Garder le tracé précédent

% Boucle pour tracer les différentes méthodes
for i = 1:length(methods)
    options.method = methods{i};  % Définir la méthode courante
    x0 = [40, 15]';  % Point initial
    
    % Appeler la fonction d'optimisation pour la méthode actuelle
    [xh, result, xval] = optimdescent(params.fonction, params, options, x0);
    
    % Tracer l'évolution du critère ou grad en fonction du nombre d'itérations
    %plot(1:length(result.crit), result.crit, [colors{i} 'o-'], 'DisplayName', methods{i});  
    plot(1:length(result.grad), result.grad, [colors{i} 'o-'], 'DisplayName', methods{i});  
    % Stocker le nombre d'itérations et le temps d'exécution
    num_iterations(i) = result.iter;
    execution_time(i) = result.time;

    % Affichage des valeurs pour déboguer
    %fprintf('Méthode: %s, Itérations: %d, Temps: %.4f\n', methods{i}, result.iter, result.time);
end

% Ajouter la légende une seule fois après avoir tracé toutes les courbes
legend(methods, 'Location', 'northeast');  % Définir la position de la légende

% Ajouter les labels et le titre
xlabel('Itérations');
ylabel('Valeur du critère');
title('Évolution de valeur du critére pour différentes méthodes');

% Activer la grille pour une meilleure lisibilité
grid on;

hold off;  % Libérer le maintien des tracés

% Afficher les résultats 
fprintf('\nRésultats de loptimisation:\n');
for i = 1:length(methods)
    fprintf('%s : Nombre ditérations = %d, Temps dexécution = %.10f secondes\n', ...
        methods{i}, num_iterations(i), execution_time(i));
end

