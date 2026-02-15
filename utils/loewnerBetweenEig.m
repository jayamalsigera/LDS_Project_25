%% Check A ⪯ B ⪯ C in the Loewner (PSD) order via eigenvalues.
%
% tf is true if (B-A) and (C-B) are PSD within tolerance tol.
%
function tf = loewnerBetweenEig(A, B, C, tol)
    if nargin < 4 || isempty(tol)
        tol = 1e-10;
    end

    BA = B - A;
    CB = C - B;

    % Symmetrize (Loewner order is for symmetric/Hermitian matrices)
    BA = (BA + BA')/2;
    CB = (CB + CB')/2;

    minEig_BA = min(eig(BA));
    minEig_CB = min(eig(CB));

    tf = (minEig_BA >= -tol) && (minEig_CB >= -tol);
end
