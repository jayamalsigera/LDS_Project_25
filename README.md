# LDS_Project_25
Robust event-triggered estimation in a sensor network. 

Version: Victor_140126

The present code simulates a sensore network, where each sensor runs an information form Kalman filter. The plant is standard statespace model of an  underdamped standard second order system. The criterion for transmitting is stochastic, and all sensors communicate with all sensors.


The code is structured in "Main.m" which contains initialization of variables, main loop and most functions. 
"StateSpacePlant.m" contains the class definition for a classic state-space model. 

Important variables are:
- Y: simlength x 1 cell array. Each cell contains a n_nodes x 1 column vector with the output as measured by each node at each time instance. Indexing is Y(time instance)
- X: simlength x 1 cell array. Each cell contains a n_states x 1 column vector with the states for each time instance. Indexing is X(time instance)
- X_hat: simlength x n_nodes cell array: Each cell contains a n_states x 1 column vector with the states as estimated by each sensor for each time instance. Indexing is X_hat(time instance, ID of the estimating node)
- q: simlength x n_nodes x n_nodes cell array: Each cell contains a nodes best knowledge of the information vector of some node for some time instance. Indexing is q(time instance,ID of node to be estimated,ID estimating node), any case where ID of node to be estimated =ID estimating node contains a nodes knowledge of itself, which is always up to date.
- Omega: simlength x n_nodes x n_nodes cell array: Each cell contains a nodes best knowledge of the inverse of the covariance matrix P used in Kalman filters for some node for some time instance. Indexing is the same as q.
- c: n_nodes x 1 array: each element is used with either 1 or 0 to denote whether a sensor decides to transmit its data to the other nodes.

Initialization consists of determining:
- simlength: How many timesteps, the simulation will run for
- n_nodes: Number of sensors
- n_states: Number of states in the statespace model
- A, B, C, D: Parameters of true statespace model. 
- A_hat, C_hat, Q, R: Parameters of Kalman filter (A_hat and C_hat are presently copies of A and C, while R is a column vector with each element i being $D(i)^2$)
- x0: Initial condition of the true statespace model.
- x0_hat: Initial state estimate of the Kalman filter (presently just a copy of x0)
- pi: n_nodes x n_nodes matrix showing the connection of each node to other nodes and with a scaling factor for fusing data
- delta: uncertainty factor used to limit the influence of bad estimates of the information vector of other nodes, when no data is received.
- epsilon: Tolerance on error between a priori state estimates and corrected state estimates.

- The variable Omega q presently initialized as $I_{n_nodes}$ for all cells Omega(1,i,j) $i,j = 1..n_nodes$
- The variable q is presently initialized as $Omega{1,i,i}\x0_hat$ for all cases where ID of node to be estimated =ID estimating node and t = 1, and a 0 n_states x 1 vector for all other cases t = 0.


