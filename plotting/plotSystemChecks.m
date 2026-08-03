function plotSystemChecks(A, B, C, sensorCount, G, tol)
% Plot PBH stabilizability, PBH detectability, and per-sensor local
% observability in a single figure.
%
% Layout:
%   Top (full width) — heatmap: which state dimensions each sensor can
%                       individually observe, via the null space of obsv(A,C_i)
%   Bottom-left      — stabilizability PBH rank at each UNIQUE unstable
%                       eigenvalue; must equal n for each λ
%   Bottom-right     — detectability PBH rank at each UNIQUE unstable
%                       eigenvalue; must equal n for each λ
%
% Parameters:
%   A           - n x n state transition matrix
%   B           - n x m input matrix
%   C           - p x n stacked output matrix (all sensors, row-ordered)
%   sensorCount - number of sensor nodes
%   G           - network digraph (node table must have isSensor column)
%   tol         - (optional) rank / null-space tolerance, default 1e-9

  if nargin < 6
    tol = 1e-9;
  end

  n = size(A, 1);

  green = [0.18 0.63 0.18];
  red   = [0.80 0.15 0.15];

  %% Deduplicate unstable eigenvalues
  eigs_A   = eig(A);
  unstable = eigs_A(abs(eigs_A) >= 1 - tol);
  [uEigs, mults] = deduplicateEigs(unstable, tol);
  nU = numel(uEigs);

  %% PBH ranks at each unique unstable eigenvalue
  stabRanks = zeros(nU, 1);
  detRanks  = zeros(nU, 1);
  for k = 1:nU
    lam = uEigs(k);
    stabRanks(k) = rank([lam*eye(n) - A, B],  tol);
    detRanks(k)  = rank([lam*eye(n) - A; C],  tol);
  end

  %% Observability heatmap: obs_map(i,j) = true if state dim j is
  %  in the observable subspace of sensor i
  obs_map = false(sensorCount, n);
  nodeIds = zeros(sensorCount, 1);
  sIdx    = 0;

  % Measurement rows per sensor (3 in this model), derived from C rather than
  % hardcoded so the slicing cannot drift out of step with the plant.
  senBlock = size(C, 1) / sensorCount;
  assert(senBlock == round(senBlock), ...
         'plotSystemChecks: size(C,1)=%d is not divisible by sensorCount=%d', ...
         size(C, 1), sensorCount);

  for i = 1:numnodes(G)
    if ~G.Nodes(i, :).isSensor; continue; end
    sIdx = sIdx + 1;
    nodeIds(sIdx) = i;
    rows = (senBlock*(sIdx-1)+1):(senBlock*sIdx);
    C_i  = C(rows, :);
    O    = obsv(A, C_i);
    N    = null(O);                   % columns span the unobservable subspace
    for j = 1:n
      ej = zeros(n, 1); ej(j) = 1;
      if isempty(N)                   % full rank → all dims observable
        obs_map(sIdx, j) = true;
      else
        obs_map(sIdx, j) = norm(N' * ej) < sqrt(tol);
      end
    end
  end
  nFullyObs = sum(all(obs_map, 2));

  %% Figure
  figure

  % --- Top: observability heatmap ---
  subplot(2, 2, [1 2])
  imagesc(double(obs_map))
  colormap(gca, [red; green])          % 0 → red, 1 → green
  clim([0 1])
  title(sprintf('Local Observability per State Dimension   (%d / %d sensors fully observable)', ...
        nFullyObs, sensorCount))
  xlabel('State dimension')
  ylabel('Sensor node')
  xticks(1:n)
  xticklabels(arrayfun(@(j) sprintf('x_{%d}', j), 1:n, 'UniformOutput', false))
  yticks(1:sensorCount)
  yticklabels(arrayfun(@(id) sprintf('S%d', id), nodeIds', 'UniformOutput', false))
  set(gca, 'TickLength', [0 0])

  % --- Bottom-left: Stabilizability ---
  subplot(2, 2, 3)
  cols = rankColors(stabRanks, n, green, red);
  b1 = bar(stabRanks, 'FaceColor', 'flat');
  b1.CData = cols;
  hold on
  yline(n, '--k', sprintf('n = %d', n), 'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left');
  hold off
  title(passFailTitle(all(stabRanks >= n), 'Stabilizability  (A, B)'))
  ylabel('rank([{\lambda}I-A, B])')
  xlabel('Unstable eigenvalue')
  xticks(1:nU)
  xticklabels(eigLabels(uEigs, mults))
  ylim([0, n + 1])
  grid on

  % --- Bottom-right: Detectability ---
  subplot(2, 2, 4)
  cols = rankColors(detRanks, n, green, red);
  b2 = bar(detRanks, 'FaceColor', 'flat');
  b2.CData = cols;
  hold on
  yline(n, '--k', sprintf('n = %d', n), 'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left');
  hold off
  title(passFailTitle(all(detRanks >= n), 'Detectability  (A, C)'))
  ylabel('rank([{\lambda}I-A; C])')
  xlabel('Unstable eigenvalue')
  xticks(1:nU)
  xticklabels(eigLabels(uEigs, mults))
  ylim([0, n + 1])
  grid on
end

%% ---- helpers -------------------------------------------------------

function [uEigs, mults] = deduplicateEigs(eigs_in, tol)
  remaining = eigs_in;
  uEigs = [];
  mults = [];
  while ~isempty(remaining)
    e    = remaining(1);
    same = abs(remaining - e) < tol;
    uEigs(end+1) = mean(remaining(same));  %#ok<AGROW>
    mults(end+1) = sum(same);              %#ok<AGROW>
    remaining = remaining(~same);
  end
end

function lbl = eigLabels(uEigs, mults)
  lbl = cell(numel(uEigs), 1);
  for k = 1:numel(uEigs)
    e = uEigs(k);
    if abs(imag(e)) < 1e-10
      es = sprintf('%.4g', real(e));
    else
      es = sprintf('%.4g%+.4gi', real(e), imag(e));
    end
    if mults(k) > 1
      lbl{k} = sprintf('\\lambda=%s (×%d)', es, mults(k));
    else
      lbl{k} = sprintf('\\lambda=%s', es);
    end
  end
end

function colors = rankColors(ranks, n, pass, fail)
  colors = zeros(numel(ranks), 3);
  colors(ranks >= n, :) = repmat(pass, sum(ranks >= n), 1);
  colors(ranks <  n, :) = repmat(fail, sum(ranks <  n), 1);
end

function t = passFailTitle(tf, label)
  if tf
    t = [label '   \color[rgb]{0.18,0.63,0.18}✓'];
  else
    t = [label '   \color[rgb]{0.80,0.15,0.15}✗'];
  end
end
