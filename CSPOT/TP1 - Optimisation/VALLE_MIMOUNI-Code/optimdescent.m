function [xh,result,xval] = optimdescent(critfun,params,options,x0)
%
% critfun : nom du fichier .m ?valuant la fonction objectif, son gradient et hessien
% params : param?tres pour l'?valuation de la fonction objectif
      %params.fonction = 'rosenbrock'
      %params.forme : fome de la fonction 'normal' ou 'moindres carrés'
% options : options n?cessaires pour la mise en oeuvre des algorithmes
%    options.method : 'gradient', 'gradient conjuge', 'Newton', 'Quasi-Newton', ...
%    options.pas : 'fixe', 'variable'
%    options.const : constante d'armijo
%    options.beta : taux de rebroussement
%    options.tolX, options.tolF, options.tolG, options.maxiter
% x0 : point initial
% xh : point final
% result : structure contenant les r?sultats
%    result.iter : nombre d'iterations r?alis?es
%    result.crit : valeur du crit?re par iteration
%    result.grad : norm du gradient par iteration 
%    result.temp : temps de calcul
%    result.stop : condition d'arret ('TolX', 'TolG', 'TolF', 'Maxiter')
% xval : valeurs des it?r?es
%


%initialization des variables
Tolx = 1;
Tolf=1;

iter = 0;

%Appel de la fonction
handle_f = str2func(critfun);
[f0, g0, H0, J0] = feval(handle_f, x0, params);
alpha_0 = 1;

crit = [f0];
grad = [norm(g0)];
xv = [];
tf = [Tolf];

f = f0;
x = x0;
g = g0;
H = H0;
d = -g0;
J = J0;

s = norm(g)^2;
B = eye(length(x0));
alpha = alpha_0;
lamb = 10^-3;
v = 10;

while (norm(g) > options.tolG) && (iter < options.maxiter) && (Tolx>options.tolX) && (Tolf>options.tolF)
tic   
        %recherche du pas optimal à chaque iter
    switch options.pas
        case 'variable'
            condition = 10000;
            alpha = alpha_0;

            %condition d'Armijo
            while condition >= options.const*alpha*g'*d
                alpha = options.beta*alpha;
                x_b = x + alpha*d;
                [f_b, g_b, H_b] = feval(handle_f, x_b, params);
                condition = f_b - f;
            end
        case 'fixe'
            alpha = alpha_0;
        otherwise
            warning('Unexpected option de pas')
   end

    switch options.method
        case 'gradient'

            x_suiv = x + alpha*d;
            [f_suiv, g_suiv, H_suiv] = feval(handle_f, x_suiv, params);
            

            d = -g_suiv;

        case 'gradient conjuge'
            x_suiv = x + alpha*d;
            [f_suiv, g_suiv, H_suiv] = feval(handle_f, x_suiv, params);
            t = s;
            s = norm(g_suiv)^2;


            d = -g_suiv + (s/t)*d;

           

        case 'Newton'
             
             if strcmp(options.pas, 'fixe')
                 alpha = 1;
             end

             x_suiv = x + alpha*d;
             [f_suiv, g_suiv, H_suiv] = feval(handle_f, x_suiv, params);
             d = linsolve(H, -g);


        case'Quasi-Newton'

             d = -B*g;

             d = alpha*d;
            
             x_suiv = x + d;
             [f_suiv, g_suiv, H_suiv] = feval(handle_f, x_suiv, params);
             Y = - g;
             Y = Y + g_suiv;
             
             B = B + (1+ (Y'*B*Y)/(Y'*d)) / (Y'*d)*(d*d') - ((B*Y*d'+ d*Y'*B)/(Y'*d));
             
            
        case 'Gauss-Newton'
            R = J'*J;
            d = linsolve(R, -g);
            x_suiv = x + alpha*d;
            [f_suiv, g_suiv, H_suiv, J_suiv] = feval(handle_f, x_suiv, params);
            J = J_suiv;

        case 'Levenberg-Marquardt'
            R = J'*J + lamb*B;
            d = linsolve(R, -g);
            x_suiv = x + alpha*d;
            [f_suiv, g_suiv, H_suiv, J_suiv] = feval(handle_f, x_suiv, params);
            J = J_suiv;
            if f_suiv < f %Si l'ité
% ration améliore la fonction objectif,On réduit lamb pour se rapprocher de la méthode de Gauss-Newton, ce qui accélère la convergence.
                lamb = lamb/v;
            elseif f_suiv >= f %sinon, On augmente lamb pour se rapprocher de la descente de gradient, assurant ainsi une meilleure stabilité
                lamb = lamb*v;
            end

    

        otherwise
            warning('Unexpected option de')
            
        
    end
    % Prochaine itération
            if (g_suiv'*d>0) %si d n'est pas direction de descente
                d = -g_suiv;
            end

            Tolx = norm(x_suiv - x);
            Tolf = abs(f_suiv - f);

            f = f_suiv;
            x = x_suiv;
            g = g_suiv;
            H = H_suiv;
            
    
            crit = [crit, f];
            grad = [grad, norm(g)];
            xv = [xv, x];
            tf = [tf, Tolf ];
            
            iter = iter + 1;
time = toc;
end

xh = x;
result.iter = iter;
result.crit = crit;
result.grad = grad;
result.tf = tf;

if iter >= options.maxiter
    result.stop = 'Nombre de iter max atteint';
elseif norm(g0) <= options.tolG
    result.stop = 'Convergence atteinte sur la norme du gradient';
elseif Tolx <= options.tolX
    result.stop = 'Convergence atteinte sur la variation des variables';
elseif Tolf <= options.tolF
   result.stop = 'Convergence atteinte sur la variation de la fonction';
else
   result.stop = 'Aucune des conditions n a été rempli';
end
result.time = time;
xval = xv;

end


