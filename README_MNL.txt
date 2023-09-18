The instances used for the experiments are in the folder MNL_instances/. See the README.txt file in that folder for info about the instances format.

This code was tested on mac.
For running this code, you need to have CPLEX installed.

In order to run this code in CLion, you can set the CPLEX_ROOT_DIR in the CMakeLists.txt file. Here,
replace the three dots with the path to
<Install_dir>/CPLEX_Studioxxx/

After this, you should be able to run the code.

The following 3 arguments need to be provided to the script:

1. (string) path to the Input instance
2. (bool) maxC:
    - true -> solve with max capacity = n/2 (with n number of products)
    - false-> solve with max capacity = n
3. (bool) milp
    - true -> solve the MILP formulation
    - false-> use LLRS procedure from the paper.
4, (Optional) id string for the experiments

The code will print some statistics (some of which reported in the paper) in the output files

    - output_mip.csv -> for each experiment, print the following statistics
            id_string : id of the experiment
            filename  : instance used for the experiment
            bool_maxC : 1 if the cardinality-constrained version of the problem was solved, 0 otherwise (UNconstrained problem)
            bool_milp  : 1 if the MILP formulation was solved, 0 otherwise (always 1 for these experiments)
            objval    : best primal bound found within time limit
            ex_code   : cplex optimation exit status
            dual_bound: best dual bound found within time limit
            gap       : primal dual gap, i.r., (bestobjective - bestinteger) / (1e-10 + |bestinteger|)
            cpu_time
            N_nodes   : number of nodes solved in the branch-and-bound tree


    - output_LLRS.csv -> for each experiment, print the following statistics
            id_string : id of the experiment
            filename  : instance used for the experiment
            bool_maxC : whether the UNcontrained or contrained version of the problem was solved
            bool_milp  : 1 if ht MILP was solved, 0 otherwise (always 0 for these experiments)
            objval    : optimal objval
            cpu time  
            nintervals_tot : total number of parametric bounds (i.e., KP relaxations) that had to be computed
            useReduction1  : whether the first  items reduction criterion was used
            useReduction2  : whether the second items reduction criterion was used
            nb_ruled_out_tot : number of ruled out items
            closethegap_time : cpu time to solve MILP+
            blb_before_gap   : primal gap at the end of the bounding procedure
            bub_before_gap   : dual gap at the end of the bounding procedure
            bounding_proc_time : cpu time to perform the bounding procedure



