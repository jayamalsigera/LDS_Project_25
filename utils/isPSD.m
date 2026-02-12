%% Check whether a matrix M is Positive Semi Definite
function tf = isPSD(M, tol)
    if nargin < 2, tol = 1e-10; end

    % enforce symmetry (helps with numerical noise)
    M = (M + M')/2;

    % PSD test via eigenvalues
    lam = eig(M);
    tf = all(lam >= -tol);
end
