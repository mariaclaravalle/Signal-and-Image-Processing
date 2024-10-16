function [f,g,H, J] = rosenbrock(x,params)

%
% x : valeur de la variable d'optimisation
% params : structure contenant les parametres necessaires pour evaluer la fonction objectif
% f,g,h : valeur de la fonction objectif, son gradient et son hessien
%

% A COMPLETER - cas  n = 2 et b = 2

if strcmp(params.forme,'normal') && (params.n == 2) %
    f = params.b*(x(2) - x(1)^2)^2 + (1-x(1))^2;
    g = [-2*(1-x(1))-4*params.b*x(1)*(x(2)-x(1)^2); 2*params.b*(x(2)-x(1)^2)];
    H = [2-4*params.b*((x(2)-3*x(1)^2)) -4*params.b*x(1); -4*params.b*x(1) 2*params.b];
    J = [];
end

if strcmp(params.forme,'moindres carrés') && (params.n == 2)
    r1 = 1-x(1);
    r2 = sqrt(params.b)*(x(2)-x(1)^2);
    f = r1^2 + r2^2;
    J = [-1 0; -2*sqrt(params.b)*x(1) sqrt(params.b)];
    g = [-2*(1-x(1))-4*params.b*x(1)*(x(2)-x(1)^2); 2*params.b*(x(2)-x(1)^2)];
    H = [2-4*params.b*((x(2)-3*x(1)^2)) -4*params.b*x(1); -4*params.b*x(1) 2*params.b];

end


