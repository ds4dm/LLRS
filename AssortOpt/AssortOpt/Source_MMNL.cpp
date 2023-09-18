#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cplex.h"
#include "cpxconst.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
# include <assert.h>

#define NB   1005   // Maximum number of products
#define MB   105    // Maximum number of classes of consumers
#define MAXP 2000   // Maximum profit

#define MAX(a, b) (((a)>(b))?(a):(b))      // max function
#define MIN(a, b) (((a)<(b))?(a):(b))      // min function
#define ABS(a, b) (((a)<(b))?(b-a):(a-b))  // abs function

#define DMAXI 999999.0 // bigM (for doubles)
#define IMAXI 19999999   // bigM (for integers)

#define sprintf_s snprintf
int n;            // Number of products
double c[NB + 1]; // Product fixed cost
int vi[NB + 1]; // (int) Product preference weights (v[0]: no-purchase option)
int maxC; // maximum cardinality assortment
/*
 *  Mine
 */
double Translate[NB+1];
bool Cons_Real_idxs[NB + 1];
int nb_ruledout_lastiter = 0;
//real_idxs: reduced idx -> real idx; Real_idxs_map:  real_idx -> reduced idx
int Real_idxs[NB+1], Real_idxs_map[NB+1];
FILE *log_out;  // Print logs
FILE *csv_out;
char filename[10000]; // Input filename
int nodelimit;
bool solveLP;
bool missing_ctrs;
bool useReduction1=false;
bool useReduction2=false;
bool useNewBounds=false;
bool useRelaxBound=false;
bool useLPbound=false;
double alpha_param=0.85, alpha_init=2.;
int ncall_tot = 0; // Total number of quick sort calls.
int ncall_unique_sort = 0; // Total number of quick sort calls.int ncall_tot = 0;
int nintervals_tot=0;// Total number of times i oculd naively call quicksort.
int nb_equal_sols = 0;
int nb_worst_sols = 0; //should stay at zero
int nb_better_sols=0;
int nb_ruled_out_tot=0;
bool bool_maxC;
int max_iters = -IMAXI;
bool first_iter=false;
int first_iter_cnt, not_first_cnt;
int counts_per_lambda=-1;
int maxiterlast_mine=-1;
int heu_calls = -1;
int ex_code_MILP_MMNL = -1;
int incumbent_updated=0;
/**************/

int m; // Number of classes of consumers
double pjm[NB + 1][MB + 1]; // Profit margin product j if sold to class m
double bm[MB + 1]; // Probability consumer belogs to class m
double vjmd[NB + 1][MB + 1]; // (double) Preference weight product j for class m (v[0][m]: no-purchase option)
int vjmi[NB + 1][MB + 1]; // (int) Preference weight product j for class m (v[0][m]: no-purchase option)

double old_p_p[NB+1], old_p_c[NB+1], old_p_vd[NB+1];

double cutoffvalue;
//mine, to save best sol
//int bestsol[NB+1], nbitems_bestsol; // Only to keep track of best solution

int mode;

int nnogood;
int nbend;

int ncall = 0;
double total_blb_comp3=-1, myub_comp3=-1;
bool cons[NB + 1];
bool stampa;

int call_p_n;
double call_p_p[NB + 1], call_p_c[NB + 1], call_p_vd[NB + 1];
bool call_taken[NB + 1];

int nintervals_firstiter;
double time_boundig_pric, time_closethegap;
double gap_rel_improv;
char id_string[50] = "no-id-given"; // to identify experiment-log files
double MMNL_MILP_end_gap=-1;
int nb_fails = -1;
int stopping_reason = -1;
int start_milp_criterion = 0;

double
parametricUB(int n, double *c, double *p, double *vd, double fstep, double lstep, int incr, bool *takenb,
             bool opt, int relgap_stop, int k_cust, int iter_sub);

int readInput_mmnl(bool prt) {
    FILE *in;
    int idummy;
    in = fopen(filename, "r");
    if (in==NULL){
        printf(" *** Warning: cannot open input file %s \n", filename);
        return 1;
    }
//    prt=false;
    fscanf(in, "%i %i", &n, &m);
    if (prt) printf("\n Nr. Items   %4i \n\n", n);
    if (prt) printf("\n Nr. Classes %4i \n\n", m);

    if (prt) printf("--i-- -----v----- ---p--- ---c--- \n");

    for (int i = 0; i <= n; ++i) {
        fscanf(in, "%i ", &idummy);
        //printf("%4i ", idummy);
        for (int k=1; k<=m; ++k){
            if (i>0) fscanf(in, "%lf %lf", &vjmd[i][k], &pjm[i][k]);
            else     fscanf(in, "%lf %lf", &vjmd[i][k], &bm[k]);
            vjmi[i][k] = int(vjmd[i][k] * 1000000 + 0.5);
            if (prt) {
                printf("%11.8lf %7.2lf\n", vjmd[i][k], pjm[i][k]);
            }
        }
        if (idummy==0) printf(""); //printf("\n");
        else{
            fscanf(in, "%lf", &c[i]);
            //printf("%7.2lf\n", c[i]);
        }
    }
    fclose(in);

    printf(" Read input file \n");
    return 0;
}

// Check Cplex output code
int checkStatus(CPXENVptr env, int status) {
    char errmsg[1024];

    if (status == 0) return 0;
//    printf("Error cases: %d\n", status);
//    printf("\n CPLEX error \n");
//    CPXgeterrorstring(env, status, errmsg);
//    printf(" %s \n", errmsg);
    return status;
}

// Add coefficient Cplex constraint
inline void addCoefConstr(int *matind, double *matval, int *nzcnt, int ind, double coef) {
    if (ind == IMAXI) return;
    matind[*nzcnt] = ind;
    matval[*nzcnt] = coef;
    *nzcnt = *nzcnt + 1;
}

bool start_milp_plus(double bestub, double myub, double bestlb, int iter, double rho){
    /****
     * 0:  iter >=100
     * 1:  bestub - blb <= 0.01 #if this does not happen, likely we have bad multipliers and we won't converge
     */
    if (start_milp_criterion==0){
        if (bestub - myub < 0.1 & iter>=100) return true;
        else return false;
    }
    else if  (start_milp_criterion==1){
        if (bestub - bestlb < 0.01) return true;
        else return false;
    }
}

double clcProfitSol3(double *x, int *ind, bool prt) {
    // Compute profit of a solution defined by cplex variables x
    double ct, pt, ptt[MB + 1], vdt[MB + 1];
    int cnt;

    ct = pt = 0.0;
    cnt = 0;
    for (int k = 1; k <= m; ++k)
        vdt[k] = vjmd[0][k];
    if (prt)
        printf(" Products: \n");
    for (int j = 1; j <= n; ++j)
        if (x[ind[j]] > 0.9) {
            if (prt)
                printf(" %4i\n", j);
            ct += c[j];
            for (int k = 1; k <= m; ++k)
                vdt[k] += vjmd[j][k];
            ++cnt;
        }

    for (int k = 1; k <= m; ++k)
        ptt[k] = 0;
    for (int j = 1; j <= n; ++j)
        if (x[ind[j]] > 0.9)
            for (int k = 1; k <= m; ++k) {
                //printf(" Add %i to class %i: %lf \n", j, k, bm[k] * pjm[j][k] * vjmd[j][k]);
                ptt[k] += bm[k] * pjm[j][k] * vjmd[j][k];
            }

    pt = 0.0;
    for (int k = 1; k <= m; ++k) {
        pt += ptt[k] / vdt[k];
        //printf(" Class %i: %lf [num %lf - denom %lf] \n", k, ptt[k] / vdt[k], ptt[k], vdt[k]);
    }
    if (prt) {
        //printf("\n Sum of preferences        %10.8lf \n", vdt);
        printf(" Total costs           %10.4lf \n", ct);
        printf(" Total marginal profit %10.4lf \n", pt);
        printf(" Final profit          %10.4lf \n", pt - ct);
        printf(" Items in solution     %4i \n", cnt);
    }

    return pt - ct;
}

// Compact form of AO under MMNL by
// Bront et al. (2009): "A Column Generation Algorithm for Choice-Based Network Revenue Management"
double compForm_Bront(bool* inside, double tlim, int mcard) {
    CPXENVptr   env;
    CPXLPptr    lp;
    int         status;
    double      obj[1], rhs[1], lb[1], ub[1], matval[2 * NB + 2];
    int         matind[2 * NB + 2], matbeg[1], nzcnt;
    char        sense[1], xctype[1], name[50], *colname[1], *rowname[1];
    double      objval, x[NB * MB + 3];
    int         ind_z[NB + 1][MB + 1], ind_y[NB + 1], ind_x[MB + 1], nvar;

    printf("\n\n Compact Formulation \n\n");

    objval = 0.0;

    env = CPXopenCPLEX(&status);
    if (checkStatus(env, status)) return objval;

    lp = CPXcreateprob(env, &status, "AO");
    if (checkStatus(env, status)) return objval;

    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    CPXsetintparam(env, CPX_PARAM_THREADS, 1);
    CPXsetintparam(env, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.000000000001);
    CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.000000000001);
    CPXsetdblparam(env, CPX_PARAM_EPOPT, 0.000000001);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 0.000000001);
    CPXsetdblparam(env, CPX_PARAM_TILIM, tlim);

    CPXchgobjsen(env, lp, CPX_MAX);

    colname[0] = name;

    nvar = 0;

    // Variables z_jk = probability of j being bought by k given assortment x in [0,1]
    lb[0] = 0.0;
    xctype[0] = 'C';
    for (int j = 0; j <= n; ++j)
        for (int k = 1; k <= m; ++k) {
            if (j == 0) {
                obj[0] = 0.0;
            } else {
                obj[0] = pjm[j][k] * bm[k] * vjmd[j][k];
                if ((pjm[j][k] < 0.00000000001) | (vjmd[j][k] < 0.00000000001)) {
                    ind_z[j][k] = IMAXI;
                    continue;
                }
            }
            ub[0] = CPX_INFBOUND;

            sprintf_s(name, sizeof(char) * 12, "z%d.%d", j, k);
            status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
            if (checkStatus(env, status)) return objval;
            ind_z[j][k] = nvar++;
        }

    // assortment variables y_j
    xctype[0] = 'B';
    lb[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        obj[0] = 0.;
        ub[0] = 1.0;
        if (!inside[j]) ub[0] = 0.0;

        sprintf_s(name, sizeof(char) * 12, "y%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
        if (checkStatus(env, status)) return objval;
        ind_y[j] = nvar++;
    }

    // denominator variables x_k
    xctype[0] = 'C';
    lb[0] = 0.0;
    for (int k = 1; k <= m; ++k) {
        obj[0] = 0;
        ub[0] = CPX_INFBOUND;

        sprintf_s(name, sizeof(char) * 12, "x%d", k);
        status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
        if (checkStatus(env, status)) return objval;
        ind_x[k] = nvar++;
    }

    rowname[0] = name;
    matbeg[0] = 0;

    // ctr SUM
    sense[0] = 'E';
    rhs[0] = 1.0;
    for (int k = 1; k <= m; ++k) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_x[k], vjmd[0][k]);
        for (int j = 1; j <= n; ++j) //changed j=0 to j=1
            addCoefConstr(matind, matval, &nzcnt, ind_z[j][k], vjmd[j][k]);
        sprintf_s(name, sizeof(char) * 15, "Sum.probs%d", k);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // ctr (8) ag. 6
    sense[0] = 'L';
    rhs[0] = 1.0;
    for (int j = 1; j <= n; ++j)
        for (int k = 1; k <= m; ++k) {
            if (ind_z[j][k] == IMAXI) continue;
            nzcnt = 0;
            //xl
            addCoefConstr(matind, matval, &nzcnt, ind_x[k], vjmd[0][k]);
            //zli
            addCoefConstr(matind, matval, &nzcnt, ind_z[j][k], -vjmd[0][k]);
            //y_i
            addCoefConstr(matind, matval, &nzcnt, ind_y[j], 1);
            sprintf_s(name, sizeof(char) * 20, "(8).z%d.%dx%dy%d", j, k, k, j);
            status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
            if (checkStatus(env, status)) return objval;
        }

    // ctr z_kj <= x_k
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j)
        for (int k = 1; k <= m; ++k) {
            if (ind_z[j][k] == IMAXI) continue;
            nzcnt = 0;
            addCoefConstr(matind, matval, &nzcnt, ind_z[j][k], 1.0);
            addCoefConstr(matind, matval, &nzcnt, ind_x[k], -1.0);
            sprintf_s(name, sizeof(char) * 20, "Rel.u%d.%dx%d", j, k, j);
            status = CPXaddrows(env, lp, 0, 1, 2, rhs, sense, matbeg, matind, matval, NULL, rowname);
            if (checkStatus(env, status)) return objval;
        }

    // ctr (9) ag. 6
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j)
        for (int k = 1; k <= m; ++k) {
            if (ind_z[j][k] == IMAXI) continue;
            nzcnt = 0;
            //zli
            addCoefConstr(matind, matval, &nzcnt, ind_z[j][k], vjmd[0][k]+vjmd[j][k]);
            //y_i
            addCoefConstr(matind, matval, &nzcnt, ind_y[j], -1);
            sprintf_s(name, sizeof(char) * 20, "(8).z%d.%dx%dy%d", j, k, k, j);
            status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
            if (checkStatus(env, status)) return objval;
        }

    status = CPXmipopt(env, lp);
    if (checkStatus(env, status)) return objval;

    int ex_code= CPXgetstat(env, lp);
    ex_code_MILP_MMNL = ex_code;

    double gap;

    CPXgetmiprelgap(env, lp, &gap);
    MMNL_MILP_end_gap = gap;
    status = CPXgetobjval(env, lp, &objval);
    if (checkStatus(env, status)) return objval;
    printf("\n Optimal solution value %12.4f \n\n", objval);

    status = CPXfreeprob(env, &lp);
    if (checkStatus(env, status)) return objval;

    status = CPXcloseCPLEX(&env);
    if (checkStatus(env, status)) return objval;


    return objval;
}

double LagrRelax_with_LLRS(bool *inside, double tlim, int mcard, double fstep, double lstep) {
    CPXENVptr env;
    CPXLPptr lp;
    int status;
    double obj[1], rhs[1], lb[1], ub[1], matval[2 * NB + 2], lhs[NB + 1][MB + 1], fc[NB + 1], blb;
    int matind[2 * NB + 2], matbeg[1], nzcnt;
    char sense[1], xctype[1], name[50], *colname[1], *rowname[1];
    double objval, x[NB * MB + 3], myub, bestub;
    int ind_u[NB + 1], ind_x[NB + 1], nvar, counter[NB + 1], iterbestub, iterbestub_overall, first[NB + 1], indsub[NB + 1], indsubk;
    double lambda[NB + 1][MB + 1], lambdab[NB + 1][MB + 1], oldlambda[NB + 1][MB + 1], oldobjval[MB + 1], old_p_c[
            NB + 1][MB + 1];
    double alpha, obju[NB + 1][MB + 1];
    bool ce[NB + 1][MB + 1], keepold[MB + 1], insol[NB + 1][MB + 1], mysol[NB + 1], oldsol[
            NB + 1][MB + 1], subopt = false, oldsubopt = false;; // false; // true;
    bool heu[NB + 1];

    // c[i] = nb. of customer types having product i in consid. set
    // ce[i][k]= true if i in consid. set of cust. k
    // first[i]= smallest cust. idx considering product i
    // counter[i] = nb of customer having i in the consid. set.
    for (int i = 1; i <= n; ++i) {
        counter[i] = 0;
        first[i] = IMAXI;
        for (int k = 1; k <= m; ++k) {
            ce[i][k] = false;
            if ((pjm[i][k] > 0.00000000001) & (vjmd[i][k]>0.00000000001)) {
                ++counter[i];
                ce[i][k] = true;
                first[i] = MIN(first[i], k);
            }
        }
    }

    //  Init all the lamda(b)[i][k] lagrangian multipliers
    for (int i = 1; i <= n; ++i)
        for (int k = 1; k <= m; ++k)
            if (ce[i][k]) {
                lambda[i][k] = lambdab[i][k] = 0.0; // i.e., separate MNL models
                oldsol[i][k] = false;
                old_p_c[i][k] = DMAXI;
            }

    for (int k = 1; k <= m; ++k) {
        keepold[k] = false;
        oldobjval[k] = 0.0;
    }

    printf(" \n\n Trivial UB \n\n");

    blb = 0.0;
    bestub = IMAXI;
    iterbestub = 0;
    alpha = alpha_init;

    // Start the subgradient search: 500 iterations in total (but the limit should not get hit in general)
    double start_time = clock();
    for (int iter = 1; iter <= 500; ++iter) {
        if (((clock()-start_time)/(double)CLOCKS_PER_SEC) > 600){
            stopping_reason = 2;
            break;
        }
        myub = 0.0;

        // lhs are the gradients for the lagrangian multipliers. Init them to zero.
        for (int i = 1; i <= n; ++i)
            for (int k = 1; k <= m; ++k)
                lhs[i][k] = 0.0;

        // heu is the heurisitc solution for for the AO under MMNL: union of the optimal assortments for each customer
        for (int i = 1; i <= n; ++i)
            heu[i] = false;

        int skipped = 0;
        bool UB_improved=false;

        // Optimize each customer class separately
        int total_p_n =  0;
        nb_ruledout_lastiter = 0;
        for (int k = 1; k <= m; ++k) {
            objval = DMAXI;
            // Not init and thus fast
            double p_c[NB + 1], p_p[NB + 1], p_vd[NB + 1], ourub;
            int p_n = 0;

            // get revs, prefs and costs of cust k
            indsubk = 0;
            p_p[0] = 0.0;
            p_vd[0] = vjmd[0][k];
            p_c[0] = 0.0;
            double sum = 0;
            for (int j = 1; j <= n; ++j)
                // If item in the consid set of that customer
                if (ce[j][k]) {
                    ++indsubk; // nb items in consid set of cust. k
                    indsub[indsubk] = j;
                    ++p_n; // number of items considered by cust. k
                    p_p[p_n] = pjm[j][k] * bm[k];
                    p_vd[p_n] = vjmd[j][k];

                    // for each item in the consid set of k, add to the costs the lambda[j][kk]
                    if (k == first[j]) {
                        p_c[p_n] = c[j] / double(counter[j]);
                        for (int kk = first[j] + 1; kk <= m; ++kk)
                            if (ce[j][kk])
                                p_c[p_n] += lambda[j][kk];
                    } else {
                        p_c[p_n] = c[j] / double(counter[j]) - lambda[j][k];
                    }
                }

            // if all the costs+lambdas for a customer are the same, then we don't need to optimize
            int oldcnt = 0;
            bool same = true;
            for (int j = 1; j <= n; ++j)
                if (ce[j][k]) {
                    ++oldcnt;
                    if (ABS(old_p_c[j][k], p_c[oldcnt]) > 0.000000001) {
                        same = false;
                    }

                    old_p_c[j][k] = p_c[oldcnt];
                }
            // oldsubopt means that this customer was solved to optimality
            // at the previous iteration
            if (same & oldsubopt) {
                ++skipped;
                keepold[k] = false;
            } else {
                keepold[k] = false;
            }

            ourub = 0.0;

            // I keep the same incumbent for cust. k
            if (keepold[k]) {
                for (int i = 1; i <= p_n; ++i)
                    mysol[i] = oldsol[i][k];
                ourub = oldobjval[k];
                goto DOPO;
            }

            stampa = false;

            // Compute parametric UB for customer k given the current set of lagrtangian multipliers
            ourub = parametricUB(p_n, p_c, p_p, p_vd, fstep, lstep, 1000, mysol, subopt, 0, k, iter);
            nb_ruledout_lastiter += nb_ruled_out_tot;
            total_p_n += p_n;

            int nb_items_in_sol;
            nb_items_in_sol = 0;
            for (int jj=1; jj<=p_n; ++jj){
                if (mysol[jj])
                    nb_items_in_sol++;
            }

            // This should never happen
            if (ourub > objval + 0.0000001) {
                printf("weird ourub>objvall +... condition\n");
                printf("Iter.%i - Class %i[%lf]: objval %lf - ourub %lf \n", iter, k, bm[k], objval, ourub);
            }

            oldobjval[k] = ourub;
            for (int i = 1; i <= p_n; ++i)
                oldsol[i][k] = mysol[i];

            DOPO:;
            myub += ourub;

            // subgradients for customer k
            for (int ii = 1; ii <= indsubk; ++ii) {
                int i = indsub[ii];
                insol[i][k] = false;
                if (mysol[ii]) {
                    heu[i] = true; // Build heuristic solution to MMNL (union of products from all customers)
                    insol[i][k] = true;

                    // lhs is the gradient, computed ad in Feldman and Topal
                    if (k == first[i]) {
                        for (int kk = first[i] + 1; kk <= m; ++kk)
                            if (ce[i][kk])
                                lhs[i][kk] += 1.0;
                    } else
                        lhs[i][k] -= 1.0;
                }
            }
            continue;

            status = CPXfreeprob(env, &lp);
            if (checkStatus(env, status)) return objval;

            status = CPXcloseCPLEX(&env);
            if (checkStatus(env, status)) return objval;
        }// end loop over all k

        double xxx[NB + 1];
        int ind[NB + 1];

        ind[0] = 0;
        xxx[0] = 0.0;
        for (int i = 1; i <= n; ++i) {
            ind[i] = i;
            xxx[i] = 0.0;
            if (heu[i])
                xxx[i] = 1.0;
        }

        blb = MAX(blb, clcProfitSol3(xxx, ind, false));
        oldsubopt = subopt;
        if (myub < bestub) {
            if (start_milp_plus(bestub, myub, blb, iter, lstep)){
                subopt = true; // small improvement in UB --> solve to optimality, to close the final gap..
            }
            bestub = myub;
            iterbestub = iter;
            iterbestub_overall = iter;
            UB_improved = true;
            for (int i = 1; i <= n; ++i)
                for (int k = 1; k <= n; ++k)
                    if (ce[i][k])
                        lambdab[i][k] = lambda[i][k];
        }

        // compute the norm of the gradient
        double denom = 0.0;
        for (int i = 1; i <= n; ++i)
            for (int k = 1; k <= m; ++k)
                denom += lhs[i][k] * lhs[i][k];

        // it means that we need to start solving the parametric UB to optimality
        if (blb>myub) {
            if (subopt)
                assert (false);
            else
                subopt = true;
        }

        if (iter - iterbestub > 2) {
            //alpha *= 0.5;
            alpha *= alpha_param;
            iterbestub = iter;
        }

        printf(" Iter.%3i: current ub %10.4lf [Best LB %10.4lf .. Best UB %10.4lf] alpha %8.6lf [subopt %2i] Skipped %3i [Improved %2i]  \n",
               iter, myub, blb, bestub, alpha, subopt, skipped, UB_improved);
        printf("Bounding time: %.4lf s\t ClosethegapTime: %.4lf s\t Ruled_out: %d/%d\t denom: %lf \t myub-bestub:%lf\n", time_boundig_pric, time_closethegap, nb_ruledout_lastiter, total_p_n, denom, myub-bestub);

        max_iters = iter; //Update iteration count
        if (denom < 0.000000001) {
            // gradient has converged to zero
            if (subopt) {
                stopping_reason=0;
                printf(" denom %lf : opt %lf \n", denom, myub);
                break;
            }
            subopt = true;
        }

        // if stepsize becomes to small, then quit.
        if (alpha <=1e-8) {
            stopping_reason=1;
            subopt = true;
            break;
        }

        // fails in improving best UB
        nb_fails = MAX(iter-iterbestub_overall, nb_fails);

        // update lambda
        for (int i = 1; i <= n; ++i)
            for (int k = 1; k <= m; ++k)
                if (ce[i][k]) {
                    oldlambda[i][k] = lambda[i][k];
                    if (denom >= 0.0000001)
                        lambda[i][k] += alpha * lhs[i][k] * (bestub)/denom;
                }
    } // iter 1->500 loop

    printf("At the end: %lf \t %lf \t %lf\n", blb, cutoffvalue, myub);
    total_blb_comp3 = blb;
    return myub;
}

// Quick sort: inner routine
int QS_Partition(double *vett, int *indici, int l, int r) {
    int i, j, k;
    double pivot, t;

    pivot = vett[l];
    i = l;
    j = r + 1;

    while (true) {
        do {
            ++i;
        } while (vett[i] >= pivot && i <= r);

        do {
            --j;
        } while (vett[j] < pivot);

        if (i >= j) break;

        t = vett[i];
        vett[i] = vett[j];
        vett[j] = t;

        k = indici[i];
        indici[i] = indici[j];
        indici[j] = k;
    }

    t = vett[l];
    vett[l] = vett[j];
    vett[j] = t;

    k = indici[l];
    indici[l] = indici[j];
    indici[j] = k;

    return j;
}

// Quicksort: main call
void QS_QuickSort(double *vett, int *indici, int l, int r) {
    int j;
    ++ncall;

    if (l < r) {
        j = QS_Partition(vett, indici, l, r);
        QS_QuickSort(vett, indici, l, j - 1);
        QS_QuickSort(vett, indici, j + 1, r);
    }
}

// Solve Linear Relaxation of Knapsack Problem
double LPRelKPbz2(int p_n, double *pr, double *w, double cap, int *indici, double oldrapp, bool *taken, double *minr) {
    double tw, ub, r[NB + 1], maxr;
    int maxp, j;

    if (cap < 0){
        for (j = 1; j <= p_n; ++j)
            taken[j] = false;
        return 0;
    }

    tw = ub = 0.0;
    *minr = IMAXI;
    for (j = 1; j <= p_n; ++j) {
        taken[j] = false;
        if (w[j]==0) r[j]=IMAXI;
        else r[j] = pr[j] / w[j];
        if (r[j] >= oldrapp & pr[j]>0) {
            tw += w[j];
            ub += pr[j];
            *minr = MIN(*minr, pr[j] / w[j]);
            taken[j] = true;
        }
    }

    if (tw > cap) {
        for (j = 1; j <= p_n; ++j)
            taken[j] = false;

        QS_QuickSort(r, indici, 1, p_n);

        ub = tw = 0.0;
        for (int j = 1; j <= p_n; ++j) {
            maxp = indici[j];
            if (pr[maxp] / w[maxp] <= 0.0)
                break;

            if (w[maxp] <= cap - tw) {
                ub += pr[maxp];
                tw += w[maxp];
                taken[maxp] = true;
                *minr = pr[maxp] / w[maxp];
                if (cap - tw < 0.00000000001)
                    break;
            } else {
                ub += (cap - tw) / w[maxp] * pr[maxp];
                tw = cap;
                *minr = pr[maxp] / w[maxp];
                break;
            }
        }
    } else {
        while (tw < cap) {
            maxp = 0;
            maxr = 0.0;
            // look for max_ratio (and item) among not taken
            for (j = 1; j <= p_n; ++j)
                if (!taken[j])
                    if (maxr < r[j]) {
                        maxr = r[j];
                        maxp = j;
                    }
            // all the not taken have a negative ration --> do not insert them in KP
            if (maxp == 0) {
                break;
            }

            // insert full quantity best item
            if (w[maxp] <= cap - tw) {
                ub += pr[maxp];
                tw += w[maxp];
                taken[maxp] = true;
                *minr = pr[maxp] / w[maxp];
                if (cap - tw < 0.00000000001)
                    break;
                // the item is the critical one. Insert fractional part.
            } else {
                ub += (cap - tw) / w[maxp] * pr[maxp];
                tw = cap;
                *minr = pr[maxp] / w[maxp];
                break;
            }
        }
    }

    return ub;
}

static int CPXPUBLIC lazycallback3(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
double x[NB + 1], val;
int status;
int j, cnt;
bool lifted;

if (mode == 2) lifted = false;
else lifted = true;

status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, call_p_n - 1);
if (status != 0)
return status;

// Compute profit of integer solution
int ind[NB + 1];
for (j = 1; j <= call_p_n; ++j)
ind[j] = j - 1;

double ct, pt, vdt;

ct = pt = 0.0;
cnt = 0;
vdt = call_p_vd[0];
for (int j = 1; j <= call_p_n; ++j)
if (x[ind[j]] > 0.9) {
ct += call_p_c[j];
vdt += call_p_vd[j];
++cnt;
}

for (int j = 1; j <= call_p_n; ++j)
if (x[ind[j]] > 0.9)
pt += call_p_p[j] * call_p_vd[j];
pt /= vdt;

val = pt - ct;

// if profit better than cutoff
if (val > cutoffvalue) {
// update cutoff
cutoffvalue = val;
//printf(" updated best primal bound %lf \n", val);
// update incumbent solution
for (int j = 1; j <= call_p_n; ++j) {
call_taken[j] = false;
if (x[ind[j]] > 0.9)
call_taken[j] = true;
}
}

int cutind[NB + 1];
double cutval[NB + 1];

cnt = 0;
for (j = 1; j <= call_p_n; ++j) {
cutind[j - 1] = j - 1;
if (x[j - 1] > 0.9) {
cutval[j - 1] = 1.0;
++cnt;
} else {
if (lifted) cutval[j - 1] = 0.0;
else cutval[j - 1] = -1.0;
}
}

status = CPXcutcallbackadd(env, cbdata, wherefrom, call_p_n, cnt - 1, 'L', cutind, cutval, CPX_USECUT_FORCE);
++nnogood;

if (status != 0)
return status;

*useraction_p = CPX_CALLBACK_SET;

return 0;
}

double get_K(double t, double rho){
    // t = either tmin or tmax
    double K = log(t) / log(1+rho);
    return K;
}

double closenogood_compform_MMNL(int n, double *c, double *p, double *vd,
                                 double tlim, double pL, double pU, double p_lb, int k_cust, int iter_sub) {
    /*
     * Single class
     */
    double st = clock();
    CPXENVptr env;
    CPXLPptr lp;
    int status;
    double obj[1], rhs[1], lb[1], ub[1], matval[2 * NB + 2], blb;
    int matind[2 * NB + 2], matbeg[1], nzcnt;
    char sense[1], xctype[1], name[50], *colname[1], *rowname[1];
    double objval, x[2 * NB + 3];
    int ind_u[NB + 1], ind_x[NB + 1], nvar;

    //printf("\n\n Compact Formulation \n\n");

//    objval = 0.0;
    blb = p_lb;
    env = CPXopenCPLEX(&status);
    if (checkStatus(env, status)) return blb;

    lp = CPXcreateprob(env, &status, "AO");
    if (checkStatus(env, status)) return blb;

    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(env, CPX_PARAM_THREADS, 1);
    CPXsetintparam(env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.000000000001);
    CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.000000000001);
    CPXsetdblparam(env, CPX_PARAM_EPOPT, 0.000000001);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 0.000000001);
    CPXsetdblparam(env, CPX_PARAM_TILIM, tlim);
    if (nodelimit >0){
        CPXsetintparam(env, CPX_PARAM_NODELIM, nodelimit);
    }
    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 3);

    CPXchgobjsen(env, lp, CPX_MAX);
    status = CPXsetdblparam(env, CPX_PARAM_CUTLO, p_lb);

    colname[0] = name;

    nvar = 0;

    //  variables u, and their obj coefficients
    lb[0] = 0.0;
    xctype[0] = 'C';
    for (int j = 0; j <= n; ++j) {
        obj[0] = p[j];

        sprintf_s(name, sizeof(char) * 12, "u%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, NULL, xctype, colname);
        if (checkStatus(env, status)) return blb;
        ind_u[j] = nvar++;
    }

    // variables x, and their obj coefficients
    xctype[0] = 'B';
    lb[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        obj[0] = -c[j];
        ub[0] = 1.0;
        if (Cons_Real_idxs[j]){
            ub[0] = 1.0;
        }
        else {
            ub[0] = 0.0; // Ruled out items
        }

        sprintf_s(name, sizeof(char) * 12, "x%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
        if (checkStatus(env, status)) return blb;
        ind_x[j] = nvar++;
    }

    rowname[0] = name;
    matbeg[0] = 0;

    // Ctr (2.b) in the paper
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], vd[0]);
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], -vd[j]);
//        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        sprintf_s(name, sizeof(char) * 15, "2b.u%du0", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return blb;
    }

    // New ctr.
    sense[0] = 'L';
    rhs[0] = pU;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], 1/vd[0]);
//        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        sprintf_s(name, sizeof(char) * 15, "UB_u0_%d", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return blb;
        break;
    }

    // New ctr.
    sense[0] = 'G';
    rhs[0] = pL;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], 1/vd[0]);
//        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        sprintf_s(name, sizeof(char) * 15, "LB_u0_%d", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return blb;
        break;
    }

    // Ctr (2.c) in the paper
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], vd[0] + vd[j]);
        addCoefConstr(matind, matval, &nzcnt, ind_x[j], -vd[j]);
//        sprintf_s(name, sizeof(char) * 15, "Rel.u%dx%d", j, j);
        sprintf_s(name, sizeof(char) * 15, "2c.u%dx%d", j, j);
        status = CPXaddrows(env, lp, 0, 1, 2, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return blb;
    }

    // Ctr (2.d) in the paper
    sense[0] = 'E';
    rhs[0] = 1.0;
    nzcnt = 0;
    for (int j = 0; j <= n; ++j)
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], 1.0);
    sprintf_s(name, sizeof(char) * 15, "Sum.u");
    status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
    if (checkStatus(env, status)) return blb;

    // Missing constraints
    sense[0] = 'L';
    for (int j = 1; j <= n; ++j) {
        if (ind_u[j] == IMAXI) continue;
        rhs[0] = 1.0;
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], 1.0);
        addCoefConstr(matind, matval, &nzcnt, ind_x[j], 1.0);
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], - double(vd[0]) / double(vd[j]));
        sprintf_s(name, sizeof(char) * 20, "MissingCtrs.u%dx%d", j, j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // Cardinality constraint
    if (maxC < n) {
        sense[0] = 'L';
        rhs[0] = maxC;
        nzcnt = 0;
        for (int j = 1; j <= n; ++j)
            addCoefConstr(matind, matval, &nzcnt, ind_x[j], 1.0);
        sprintf_s(name, sizeof(char) * 15, "Card");
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return blb;
    }

    if (checkStatus(env, status)) return blb;

    // Solve Linear Relaxation
    if (solveLP){
        CPXchgprobtype(env, lp, 0);
        status = CPXlpopt(env,lp);
        // Opt Val
        status = CPXgetobjval(env, lp, &objval);
        double et = clock();
        printf("Objective Value: %f\n", objval);
        return blb;
    }

    // Solve MIP
    status = CPXmipopt(env, lp);

    if (checkStatus(env, status)) return blb;

    // Get Opt Sol Stats:
    // Termination code (optimal vs. time limit vs optimal within tolerance)
    int ex_code= CPXgetstat(env, lp);

    double gap;
    CPXgetmiprelgap(env, lp, &gap);

    // Opt Obj Val
    status = CPXgetobjval(env, lp, &objval);
    if (checkStatus(env, status)) return blb;
    double et = clock(); //mine

    double dual_bound;
    status = CPXgetbestobjval(env, lp, &dual_bound);

    int N_nodes =  CPXgetnodecnt(env, lp);

//    printf("\n Optimal solution value %12.8f \n\n", objval);

    // Get Optimal solution
    status = CPXgetx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1);
//    status = CPXgetx(env, lp, x, 0, n);
    if (checkStatus(env, status)) return blb;


    if (objval > blb){
        // update incumbent solution
        // Compute profit of integer solution
        incumbent_updated = 1;
        int ind[NB + 1];
        for (int j = 1; j <= n; ++j)
            ind[j] = j;// - 1; // Since in x i don't have the zero

        for (int j = 1; j <= n; ++j) {
            call_taken[j] = false;
            if (x[ind_u[j]] > 0.000000001)
                call_taken[j] = true;
        }
    }
//    else printf("Not updating: objval= %.5lf\t blb = %.5lf\n", objval, blb);

    status = CPXfreeprob(env, &lp);
    if (checkStatus(env, status)) return blb;

    status = CPXcloseCPLEX(&env);
    if (checkStatus(env, status)) return blb;

    return MAX(blb, objval); //objval;
}

double
parametricUB(int n, double *c, double *p, double *vd, double fstep, double lstep, int incr, bool *takenb,
             bool opt, int relgap_stop, int k_cust, int iter_sub){
    /*
     * fstep: first level of grid density (granularity rho in the paper: how many digits)
     * lstep: last  level of grid density
     * incr = 50: nb of intervals used to split [pmin, pmax] after the last iteration, for the "enumerative phase"
     */
    double lbs, *ubs, bestlb, bestub, *rappmin;
    int denomin, denomout, indici[NB+1], maxiterlast;
    double pmax, pmin;
    double profitti[NB + 1], pesi[NB + 1];
    // taken is the vector used to memorize current best lb/incumbent
    bool taken[NB + 1], wcard; //wcard: wether cardinality ctr is not satisfied.

    double st, et, bounding_start;
    double rho = fstep;
    bool firstiter = true;
    int old_p_n ;
    st = clock();
    bounding_start = clock();

    bestlb = 0.0;
    wcard = true;
    nb_ruled_out_tot = 0;
    incumbent_updated = 0;

    // To remember the original item idx, after reductions (and compare solutions with other methods)
    for (int i=0; i<=n; ++i) {
        Real_idxs[i]=i;
        Real_idxs_map[i]=i;
        if (i!=0) {
            Cons_Real_idxs[i] = true;
            old_p_c[i] = c[i];
            old_p_p[i] = p[i];
            old_p_vd[i] = vd[i];
        }
    }
    old_p_n = n;

    // Aux vector to store solution for this customer
    for (int i = 1; i <= n; ++i) {
        takenb[i] = false;
    }

    bool zero = true;
    // Naive solution
    for (int i = 1; i <= n; ++i)
        if (p[i] * vd[i] / (vd[i] + vd[0]) > c[i]) {
            zero = false;
            break;
        }

    if (zero)
        return bestlb;

    // Compute pmin, pmax
    double pref_sum = 0;
    double pref_min = IMAXI;
    for (int i =1; i <= n; ++i) {
        pref_sum += vd[i];
        if (vd[i] > 0) pref_min = MIN(pref_min, vd[i]);
    }

    // Starting range, as in Feldman and Topaloglu (2015) and Kunnumkal (2019) (Unnormalied v0, as in the latter)
    pmin = 1/(vd[0] + pref_sum);
    pmax = 1/(vd[0] + pref_min);

    double last_gap = 100;
    double last_ub = IMAXI;
    AGAIN:; // for every rho in [fstep, lstep]:

    // Get exponents K
    int k_min = int(get_K(pmin, rho));
    if (k_min<0) k_min--;
    int k_max = int(get_K(pmax, rho));
    if (k_max>0) k_max++;

    if (pow(1+rho, k_min) > pmin){
        // Should never happen
        printf("bestlb: %lf \t best ub: %lf \t fstep:%.15lf \t rho:%lf, \t k_cust: %d\n", bestlb, bestub, fstep, rho, k_cust);
        printf("pmin: %.15lf\t pow:%.15lf\t kmin:%d\tkmax:%d\n", pmin, pow(1+rho, k_min), k_min, k_max);
        assert (false);
    }
    assert (pow(1+rho, k_min) <= pmin);
    assert (pow(1+rho, k_max) >= pmax - 1e-8);

    if (firstiter) {
        nintervals_firstiter += k_max - k_min + 1;
        firstiter = false;
    }

    maxiterlast = 0;

    // Keep track over lower/upper bound for each denomitor (useful for "pruning" the [pmin,pmax] range)
    ubs = (double *) malloc((k_max -k_min + 1) * sizeof(double));

    // Save critical ratios, i.e., p_s/w_s (with s critical item)
    // for each LP-relaxed KP corresponding to solutions in [pL,pU]
    rappmin = (double *) malloc((k_max -k_min + 1) * sizeof(double));

    bestub = 0.0;

    // KP weights
    for (int j = 1; j <= n; ++j)
        pesi[j] = vd[j];

    ncall = 0; // Number of Calls to Quick sort per iteration
    int olditer, newiter; //helper vars for cardinality ctr case

    // The algorithm implementation needs increasing capacity, so we traverse the range from k_max to k_min
    for (int k = k_max; k >= k_min+1; --k) {
        double pL= pow(1+rho, k-1);
        double pU= pow(1+rho, k);

        // Unsorted items indexes ([1,2,3,...,N])
        for (int j = 1; j <= n; ++j)
            indici[j] = j;

        // Total product-specific profit for a given denom (t)
        profitti[0] = 0.0;
        for (int j = 1; j <= n; ++j) {
            profitti[j] = pU * p[j] * vd[j] - c[j];
        }

        double minr; // critical ratio p_s/w_s for the KP corresponding to a given interval (still to compute)
        double oldrapp; // critical ratio p_s/w_s for the KP corresponding to the previous interval
        if (k == k_max) oldrapp = -IMAXI;
        else oldrapp = rappmin[(k) - k_min];

        // Compute upperbounds. wcard=True if the maxiterlast==0, which means the lagrangian multiplier=0
        double newub;
        newub = LPRelKPbz2(n, profitti, pesi, 1 / pL - vd[0], indici, oldrapp, taken, &minr);

        ubs[k-1-k_min] = newub;

        rappmin[k-1-k_min] = minr;

        lbs = 0.0;
        // update
        if (ubs[k-1-k_min] > bestlb) {
            double ct, pt, vdt;

            ct = pt = 0.0;
            vdt = vd[0];
            for (int i = 1; i <= n; ++i)
                if (taken[i]) {
                    ct += c[i];
                    vdt += vd[i];
                    pt += p[i] * vd[i];
                }
            pt /= vdt;
            lbs = pt - ct;
        }

        bestub = MAX(bestub, ubs[k-1-k_min]);

        if (lbs>bestlb) {
            bestlb = MAX(bestlb, lbs);
            for (int i = 1; i <= n; ++i)
                takenb[i] = taken[i];
        }
    }
    et = clock();
    nintervals_tot += k_max - k_min + 1;

    // Prune (restrict) the range of intervals over which we compute the bounds.
    denomin = denomout = 0;
    pmin = pmax = 0.;
    int new_kmin,new_kmax;
    for (int k = k_min; k <= k_max-1; ++k) {
        double pL= pow(1+rho, k);
        double pU= pow(1+rho, k+1);
        if (ubs[k-k_min] < bestlb)
            ++denomout;
        else {
            ++denomin;
            if (pmin == 0) {
                pmin = pL; //  minin= smallest t: z(t, pho)> bestlb
                new_kmin= k;
            }
            pmax = pU; //  maxin= biggest t: z(t, pho)> bestlb
            new_kmax= k+1;//k+1
        }
    }
    k_min = new_kmin;
    k_max = new_kmax;

    // Reduction procedures
    for (int jj = 0; jj <= n; ++jj)
        cons[jj] = true;

    int cntout = 0;

    if ((!useReduction1)) goto SKIPRED1;

    // First reduction
    for (int j = 1; j <= n; ++j) {
        double prof;
        prof = pmax * p[j] * vd[j] - c[j];

        if (prof > -0.1)
            goto TIENI;

        cons[j] = false;
        Cons_Real_idxs[Real_idxs[j]]= false;
        p[j] = 0.0;
        c[j] = 1000000.0;
        ++cntout;
        ++nb_ruled_out_tot;

        TIENI:;
    }
    SKIPRED1:;
    if ((!useReduction2)) goto SKIPRED;

    // Second Reduction procedure
    cntout = 0;
    for (int j = 1; j <= n; ++j) {
        if (!cons[j]) continue;

        for (int k = k_min; k <= k_max-1; ++k) {
            double prof;
            double pL= pow(1+rho, k);
            double pU= pow(1+rho, k+1);
            prof = pU * p[j] * vd[j] - c[j];

            // Compute the bound only for products such that x_j=0 when computing the linear relax z(t, rho)
            if (prof / vd[j] >= rappmin[(k) - k_min])
                goto TIENI2;

            // Compute upper bound for x[j] fixed to 1.
            if (ubs[k-k_min] - rappmin[(k)-k_min] * vd[j] + prof > bestlb - 1) {
                goto TIENI2;
            }
        }
        ++nb_ruled_out_tot;
        Cons_Real_idxs[Real_idxs[j]] = false;
        cons[j] = false;
        p[j] = 0.0;
        c[j] = 1000000.0;
        ++cntout;

        TIENI2:;
    }
    SKIPRED:;
    int counter = 0;

    // Save old solution, and zero out new one
    for (int i = 1; i <= n; ++i) {
        call_taken[i] = takenb[i];
        takenb[i] = false;
    }

    // Update problem variables
    for (int j = 1; j <= n; ++j) {
        if (cons[j]) {
            ++counter;
            p[counter] = p[j];
            c[counter] = c[j];
            vd[counter] = vd[j];
            Real_idxs[counter] = Real_idxs[j];
            Real_idxs_map[Real_idxs[j]] = counter;
            if (call_taken[j])
                takenb[counter]= true;
        }
    }
    n = counter;
    for (int j = n + 1; j <= old_p_n; ++j)
        takenb[j]=false;
    ///////////////////////////////////////

    double dual_gap = ABS((bestub - bestlb)/bestlb+1e-10, 0) * 100;
    if (relgap_stop==1){
        // bound improv. Theoretically this should be low
        gap_rel_improv = ABS((last_ub - bestub)/last_ub*100, 0);
    }
    else{
        // Put emphasis on outputting high quality solutions
        gap_rel_improv = ABS((dual_gap - last_gap)/last_gap*100, 0);
    }

    last_gap = dual_gap;
    last_ub = bestub;

    if (((relgap_stop>0) & (gap_rel_improv>1)) | ((relgap_stop==0) & (rho > lstep+1e-8) )){
        if (maxiterlast == 0)
            wcard = false;
        else
            wcard = true;

        free(ubs);
        free(rappmin);

        rho /= 10.0;
        goto AGAIN;
    }
    // End of boudning procedure

    double bounding_end = clock();
    time_boundig_pric += (bounding_end-bounding_start)/(double)CLOCKS_PER_SEC;


    mode = 1;
    int n_after_red = n;
    n = old_p_n;
    // set all variables to original value
    for (int i=1; i <= n; ++i) {
        c[i] = old_p_c[i];
        p[i] = old_p_p[i];
        vd[i] = old_p_vd[i];
    }


    // taken_b has already been updated to be equal to taken (line 1344, or close)
    // so taken_b already has best incumbent
    // set call_taken = takenb, since call_taken will be modified in closenogood_compform_MMNL
    for (int i = 1; i <= n; ++i) {
        call_taken[i] = false;
    }
    for (int i=1; i<= n_after_red; ++i){
        call_taken[Real_idxs[i]] = takenb[i];
    }

    // mode == 1: MILP+ to close the gap
    cutoffvalue = bestlb;
    assert (mode==1);
    if (opt) {
        st = clock();
        cutoffvalue = closenogood_compform_MMNL(n, c, p, vd, 600., pmin, pmax, cutoffvalue, k_cust, iter_sub);

        time_closethegap += (clock() - st) / (double) CLOCKS_PER_SEC;
    }

    for (int i = 1; i <= n; ++i)
        takenb[i] = call_taken[i];



    free(ubs);
    free(rappmin);

    return cutoffvalue;
}


int main(int argc, char* argv[]) {
    double st, et;
    int inst_type = 0;
    bool inside[NB + 1];

    // Input filename
    strcpy(filename, argv[1]);

    double fstep, lstep;
    fstep = 0.1;
    lstep = 0.0001;
    alpha_init = 1; //Initial stepsize for lagrangian relaxation
    alpha_param = 0.5; //alpha = alpha * alpha_param when UB is not improved
    // if (mip), solve bront et al.
    bool mip = (strcmp(argv[2],"true")==0 | strcmp(argv[2],"True")==0) ? 1 : 0;
    if (argc >= 4) {
        strcpy(id_string, argv[3]); //id string of the experiment
    }
    else{
        strcpy(id_string, "No_id_string");
    }
    printf("Here 2\n");
    useReduction1 = true;
    useReduction2 = true;
    //parameter to decide when to start running MILP+ to solve AOPC for each customer type, rather than just using the
    // bounding procedure
    start_milp_criterion = 1;

    printf("Parameters:");
    printf("fstep: %.15lf \t lstep:%.15lf\n", fstep, lstep);
    printf("alpha_init: %.4lf \t alpha_param:%.4lf\n", alpha_init, alpha_param);
    printf("useRed1: %d \t useRed2:%d\n", useReduction1, useReduction2);

    // cardinality case or not: Non cardinality constraint was studied for the AOP under MMNL
    bool_maxC = false;

    // whether to only solve the linear relaxation when closing the gap
    solveLP = false;

    char csv_name[1000];
    strcpy(csv_name,"../output_MMNL.csv");

    if ((csv_out = fopen(csv_name, "a")) == NULL) {
        printf(" *** Warning: cannot open csv file %s\n", csv_name);
        return 1;
    }

    nodelimit = -1;


    if (readInput_mmnl(false)) return -1;
    printf("Nb customer types: %d\n", m);
    fprintf(csv_out, "%s, %s, %d, %d,", filename, id_string, bool_maxC, mip);

    for (int i = 0; i <= n; ++i)
        inside[i] = true;

    if (bool_maxC) maxC = int(n / 2);
    else maxC = n;

    nnogood = nbend = 0;
    mode = 3;

    st = clock();
    if (!mip) myub_comp3 = LagrRelax_with_LLRS(inside, 40.0, maxC, fstep, lstep);
    else myub_comp3 = compForm_Bront(inside, 600.0, maxC);
    cutoffvalue = total_blb_comp3; // just for printing the same as in ParametricUB
    et = clock();

    printf(" Best LB %10.4lf \n", cutoffvalue);
    printf(" CPU time %lf \n", (et - st) * 0.001);
    fprintf(csv_out, "%f, %f, %f,", cutoffvalue, myub_comp3, (et-st)/(double)CLOCKS_PER_SEC);
    // New
    fprintf(csv_out, "%d,", max_iters);
    // New
    fprintf(csv_out, "%f, %f, %d, %.10lf,%d,", alpha_init, alpha_param, ex_code_MILP_MMNL, MMNL_MILP_end_gap, stopping_reason);
    fprintf(csv_out, "%lf, %lf\n", time_closethegap, time_boundig_pric);

    FILE *out;
    //error_t err; //errno_t, error_t
    if ((out = fopen("out.txt", "a")) == NULL) {
        printf(" > Cannot open output file \n");
        exit(0);
    } else {
        fprintf(out, "%10.4lf\t%7.2lf\n", cutoffvalue, (et - st) * 0.001);
        fclose(out);
    }
    printf("maxiters: %d\n", max_iters);
    printf("Nb sorting at first iter: %d\n", first_iter_cnt);
    printf("Nb sorting NOT at first iter: %d\n", not_first_cnt);
    printf("Nb sorting total: %d\n", ncall_unique_sort);
    printf("Nb sorting per lambda: %d\n", counts_per_lambda);
    printf("Nb ruled out tot: %d\n", nb_ruled_out_tot);
    printf("Nb intervals tot: %d\n", nintervals_tot);
    printf("Nb ruled out LAST ITER: %d\n", nb_ruledout_lastiter);
    printf("EX CODE MILP MMNL: %d\n", ex_code_MILP_MMNL);
    printf("MILP GAP: %.10lf\n", MMNL_MILP_end_gap);
    printf("Stopping reason: %d\n", stopping_reason);
    printf("Nb fails: %d\n", nb_fails);
    fclose(csv_out);
}
