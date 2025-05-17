-------------------------------------- CS566 ASSIGNMENT 6------------------------------------------------
NAME: Lakshya Onkara
ROLL NO.: 210101062

AIM: Solution to HMM problem-3
EXECUTION: Build and run the code on Visual Studio 2010. Use F5 key to run the code.

INPUT: Initial_Model.txt, HMM_OBSERVATION_SEQUENCE_1.txt
OUTPUT: Optimal model parameters

CONSTANTS:
1. M: The number of observation symbols per state (32)
2. N: The number of states (5)
3. MAX_T: Maximum number of observations (160)

VARIABLES:
1. A[N+1][N+1]: Array to store the state transition probability distribution.
2. B[N+1][M+1]: Array to store the observation symbol probability distribution.
3. Pi[N+1]: Array to store the initial state distribution.
4. T: Variable to store number of observations
5. O[MAX_T+1]: Array to store the observation sequence.

6. alpha[MAX_T+1][N+1]: Array to store the forward variable values.
7. beta[MAX_T+1][N+1]: Array to store the backward variable values.
8. P_O_given_lambda: Long double variable to store the solution of problem 1.

9. delta[MAX_T+1][N+1]: Array to store the best score along a single path, at time t, which accounts for the first t observations and ends in state i.
10. psi[MAX_T+1][N+1]: Array to keep track of the argument that maximized delta[t+1][j] for each t and j.
11. q_star[MAX_T+1]: Array to store the optimal state sequence.
12. p_star: long double to store the solution of problem-2.

13. gamma[MAX_T+1][N+1]: Array to store the gamma values.
14. xi[MAX_T+1][N+1][N+1]: Array to store the probability of being in state Si at time t & Sj at time t+1 given the observation sequence & the model.
15. A_bar[N+1][N+1]: Array to store the reestimated state transition probability distribution.
16. B_bar[N+1][M+1]: Array to store the reestimated observation symbol probability distribution.
17. Pi_bar[N+1]: Array to store the reestimated initial state distribution.

FUNCTIONS:
1. computeXi():
	Function to compute the probability of being in state Si at time t & Sj at time t+1 given the observation sequence & the model.
2. computeGamma():
	Function to compute the probability of being in state Si at time t given the observation sequence & the model.
3. reestimateParameters():
	Function to iteratively update and improve the HMM parameters.

PROCEDURE:
1. Read the input files and populate the A, B, Pi and O arrays.
2. Call the functions.
3. Echo p_star sfter each iteration.