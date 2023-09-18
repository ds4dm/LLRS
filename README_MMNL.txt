The instances used for the experiments are in the folders MMNL_25, MMNL_50, MMNL_75. the number in the folder denotes
the number of customer types used to generate the corresponding set of instances.

These instance were generated according to "Bounding optimal expected revenues for assortment optimization under mixtures of multinomial logits"
by Feldmand and Topaloglu (2015). In particular, we used the same code kindly provided by the author, and converted the instances to the format
that follows.

The full data generation process and paraneters are summarized in the paper.

The following parameters were used for naming the instances:
    - G : number of customer types
    - K : maximum possible value of the parameter k_j, used to vary the the magnitude of the preference weights
    - P0 (gamma): maximum probability of the no-purchase option
    - S : number of staple products (always 40 in our experiments)
    - seed: random seed used to generate the instance

Instance name: MMNL_<G>_<K>_<P0>_<S>_<seed>.txt
   e.g.        MMNL_25_5.0_1.0_40_91.txt = Instance with 25 customer types, K=5, P0=1, 40 staple products and random seed equal to 91

let
    - v[j][k], denote the preference weight customer k for item j
    - p[j], c[j], denote the profit and cost of item j (same for all customers)
    - let b[k] denote the probability of a customer belonging to type k

Each instance contains the following informations:

n , G -> number of products (always 100 in our experiments), number of customer types
0 v[0][1] b[1] v[0][2] b[2] .... v[0][G] b[G]
1 v[1][1] p[1] v[1][2] p[1] .... v[1][G] p[1] c[1]
.
.
.
j v[j][1] p[j] v[j][2] p[1] .... v[j][G] p[j] c[j]
.
.
.
n v[n][1] p[n] v[n][2] p[n] .... v[n][G] p[n] c[n]

note the c[j] is always 0 in the instances we provide.

This code was tested on mac.
For running this code, you need to have CPLEX installed.

In order to run this code in CLion, you can set the CPLEX_ROOT_DIR in the CMakeLists.txt file. Here,
replace the three dots with the path to
<Install_dir>/CPLEX_Studioxxx/

After this, you should be able to run the code.

The following 2 arguments need to be provided to the script:

1. (string) path to the Input instance
2. (bool) milp
    - true -> solve the MILP formulation
    - false-> use LLRS procedure from the paper.
3. (Optional) id string for the experiments

The code will print some statistics (some of which reported in the paper) in the output files

    - output_MMNL.csv -> for each experiment, print the following statistics
            id_string : id of the experiment
            filename  : instance used for the experiment
            bool_maxC : 1 if the cardinality-constrained version of the problem was solved, 0 otherwise (UNconstrained problem)
            bool_milp  : 1 if the MILP formulation was solved, 0 otherwise (always 1 for these experiments)
            objval    : best primal bound found within time limit
            dual_bound: best dual bound found within time limit
            cpu_time  : total cpu time in seconds to run whether the MILP formulation or LLRS
            nb_iters   : nb of iterations of the lagrangian relaxation
            alpha_init: initial stepsize to update lagrangian multipliers (default is 1)
            alpha_param: factor to decrease the stepsize (default is 0.5)
            ex_code_MILP_MMNL : cplex optimation exit status (-1 if LLRS is executed rather than MILP formulation)
            MMNL_MILP_gap : optimality gap at time limit after executing MILP formulation (-1 if LLRS is executed)
            ex_code_LLRS : exit code of LLRS
                - 0 if optimal solution was found
                - 1 if stepsize becomes to small (< 1e-8)
                - 2 if time limit is reached
            time_milp_plus : total time (seconds) spent solving MILP+
            time_bounding_proc: total time (seconds) spent solving the bounding procedure



