%% Robust Distributed Kalman Filter
%
% Implementation based on Ghion, Zorzi (2023).
%
classdef RDKF
  properties
    Ts
    T
    % Tunning Parameters
    alpha
    beta
    delta
    b
    % Network Graph and Parameters
    G
    Pi
    N
    S
    % State estimate history
    X_hat
    % Model Matrices
    A
    C
    Q
    R
    n
    % Stats
    RMSE
    txRate
  end

  methods
    function self = RDKF(plant, Ts, T, G, alpha, beta, delta,b)
      self.Ts = Ts;
      self.T = T;
      self.G = G;

      self.alpha = alpha;
      self.beta = beta;
      self.delta = delta;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      self.Pi = calcFusionWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.b=b;

      self.X_hat = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
      % It's unrealistic to assume all nodes share same initial conditions
      % (Battistelli & Chisci, 2014), specially with "perfect knowledge", but

      q_pred = repmat(P0 \ x0_hat, 1, self.N);
      Psi = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      % Initializing the "global" predictions, assuming c_t = 1 for all nodes in the first
      % iteration (i.e. the first fusion step only relies on the local filtered values).
      q_bar = nan(self.n, self.N);
      Psi_bar = nan(self.n, self.n, self.N);
        %Iterate through all time steps
      for t = 2:self.T + 1
        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(q_pred, Psi, y); %Correction step

        for i = 1:self.N
          self.X_hat(:, i, t) = pinv(Omega_upd(:, :, i)) * q_upd(:, i); %Return to state spae form from information form
        end

        c_t = self.exchange(self.X_hat(:, :, t), Omega_upd, q_bar, Psi_bar); %Evaluate trigger condition
        self.txRate(t) = sum(c_t) / self.N; %Evaluate transmission rate, ie. what fraction of nodes transmit.



        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Psi_bar); %Fuse data
        [q_pred, Psi] = self.getLocalPriors(q_fused, Omega_fused); %Prediction steps for transmitting nodes
        [q_bar, Psi_bar] = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Psi_bar); %Prediction step for non-transmitting nodes
         %t
      end

      self.RMSE = self.calculateRSME(self.X_hat, X);
    end

    %% Correction/Update/Measurement step
    % Update the local information pair to obtain q_k|k and Omega_k|k of each
    % node
    function [q_upd, Omega_upd] = update(self, q_pred, Psi, y)
      q_upd = zeros(self.n, self.N);
      Omega_upd = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          idx = (2 * i - 1):(2 * i);

          y_i = y(idx);
          C_i = self.C(idx, :);
          R_i = self.R(idx, idx); % R (block) diagonal since D (block) diagonal

          q_upd(:, i) = q_pred(:, i) + C_i' * (R_i \ y_i);
          Omega_upd(:, :, i) = Psi(:, :, i) + C_i' * (R_i \ C_i);
        else
          q_upd(:, i) = q_pred(:, i);
          Omega_upd(:, :, i) = Psi(:, :, i);
        end
      end
    end

    %% Information Exchange
    % Evaluate whether to transmit measurement
    function c_t = exchange(self, X_hat, Omega_upd, q_bar, Psi_bar)
      c_t = ones(self.N, 1);
      if any(isnan(Psi_bar), "all")
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        Omega_i = Omega_upd(:, :, i);

        x_bar_i = Psi_bar(:, :, i) \ q_bar(:, i);

        e = X_hat(:, i) - x_bar_i; % Discrepancy from Prediction since last transmission
        eNorm = e' * Omega_i * e; % Weighted Euclidean Norm

        lower = (1 / (1 + self.beta)) * Omega_i;
        upper = (1 + self.delta) * Omega_i;

        if eNorm <= self.alpha && loewnerBetweenEig(lower, Psi_bar(:, :, i), upper) %Deterministic trigger condition
          c_t(i) = 0;
        end
      end
    end

    %% Information Pair Fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_fused = zeros(self.n, self.N);
      Omega_fused = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        [~, nids] = inedges(self.G, i);

        for j = nids'
          pi_ij = self.Pi(i, j);

          if (i == j) || c_t(j)
            % Node i has received from j or this is a self-loop (node has access to its own local info)
            q_fused(:, i) = q_fused(:, i) + pi_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_upd(:, :, j);
          else
            q_tilda = (1 / (1 + self.delta)) * q_bar(:, j);
            Omega_tilda = (1 / (1 + self.delta)) * Psi_bar(:, :, j);

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_tilda;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_tilda;
          end
        end
      end
    end

    %% Prediction Step
    %Robust prediction for transmitting nodes
    function [q_pred, Psi] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);
      Psi=zeros(self.n, self.n, self.N);
      for i = 1:self.N
        q_i_F = q_fused(:, i);
        Omega_i_F = Omega_fused(:, :, i);

        Omega_pred(:, :, i) = self.updateOmega(Omega_i_F);

         theta=self.findTheta(Omega_pred(:, :, i)); %Find theta satisfying gamma(Omega,theta)=b in Ghion, Zorzi (2023)
        Psi(:,:,i)=Omega_pred(:, :, i)-theta*eye(self.n); %Robust information set: Psi
        q_pred(:, i) = Psi(:, :, i) * self.A * (Omega_i_F \ q_i_F); %Robust information set: q
      end
    end
    %Propagation of unupdated information sets and perform robust
    %prediction
    function [q_bar_next, Psi_bar_next] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_bar_next = zeros(self.n, self.N);
      Omega_bar_next = zeros(self.n, self.n, self.N);
      Psi_bar_next=zeros(self.n, self.n, self.N);
      for i = 1:self.N
        if c_t(i)
          q_i_check = q_upd(:, i);
          Omega_i_check = Omega_upd(:, :, i);
        else
          q_i_check = q_bar(:, i);
          Omega_i_check = Psi_bar(:, :, i);
        end

        Omega_bar_next(:, :, i) = self.updateOmega(Omega_i_check);
        theta_bar=self.findTheta(Omega_bar_next(:, :, i)); %Find theta satisfying gamma(Omega,theta)=b in Ghion, Zorzi (2023)
        Psi_bar_next(:,:,i)=Omega_bar_next(:,:,i)-theta_bar*eye(self.n); %Robust information set: Psi
        q_bar_next(:, i) = Psi_bar_next(:,:,i) * self.A * (Omega_i_check \ q_i_check); %Robust information set: q
      end
    end

    function newOmega = updateOmega(self, Omega)
      invQ = self.Q \ eye(self.n);
      invQA = invQ * self.A;
      % TODO: Review variable name
      foo = (Omega + (self.A' * invQA)) \ eye(self.n);

      newOmega = invQ - invQA * foo * invQA'; % Assuming Q = Q'
    end


    %Creat and solve equation to find theta
        function theta = findTheta(self, Omega)
        % Compute eigenvalues of inv(Omega)
        lambda = eig(Omega);
        I=eye(self.n);
        invOmega=Omega\I;
        % Define y(theta)
        yfun = @(th) 0.5 * ( trace(inv(I - th*invOmega) - I) ...
               + log(det(I - th*invOmega)) );

        % Upper bound from theory
        lambda_min = min(lambda);


        % Solve y(theta) = b using bisection
        theta = self.bisect(@(th) yfun(th) - self.b, 0, 0.999*lambda_min, 1e-6);
        end

    %Bisection method used to find theta
    function theta = bisect(self,fun, a, b, tol)
    fa = fun(a);

    while (b - a) > tol
        m = 0.5*(a + b);
        fm = fun(m);

        if fa * fm <= 0
            b = m;
        else
            a = m;
            fa = fm;
        end
    end

    theta = 0.5*(a + b);
end


    %% RMSE Calculation
    % TODO: Maybe we can move this to `utils`?
    function [rmse] = calculateRSME(self, X_hat, X)
      rmse = zeros(self.T + 1, 1);
      for t = 1:self.T
        err = X_hat(:, :, t) - X(:, t); % implicit expansion over N
        rmse(t) = sqrt(mean(sum(err .^ 2, 1)));
      end
    end

    %% Plotting
    function plotTrajectory(self, X)
      % TODO: Would be cool if we could plot P(t) somehow
      % TODO: Restrict axis to ranges of X

      meanX_hat = mean(self.X_hat, 2);

      figure
      plot(meanX_hat(3, :), meanX_hat(4, :));
      hold on
      plot(X(3, :), X(4, :));
      hold off
      title("RDKF Estimated Trajectory")
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({"RDKF", "Actual Model"})
      grid()
    end
  end
end
