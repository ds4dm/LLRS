#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ilcplex/cplex.h"
#include "ilcplex/cpxconst.h"
#include <assert.h>

#define NB   1005   // Maximum number of products
#define MAXP 2000   // Maximum profit

#define MAX(a, b) (((a)>(b))?(a):(b))      // max function
#define MIN(a, b) (((a)<(b))?(a):(b))      // min function

#define IMAXI 19999999   // bigM (for integers)

#define sprintf_s snprintf

int n;            // Number of products
double p[NB + 1]; // Product profit margin
double c[NB + 1]; // Product fixed cost
double vd[NB + 1]; // (double) Product preference weights (v[0]: no-purchase option)
int vi[NB + 1]; // (int) Product preference weights (v[0]: no-purchase option)
int maxC; // Maximum cardinality of assortment

// Input/Output files
FILE *csv_out;
char filename[10000]; // Input filename
char id_string[50] = "no-id-given"; // to identify experiment-log files

// Input arguments
int nodelimit=-1; //if nodelim <0 --> full MILP (compform)
bool solveLP=false; //if true, compute the linear relaxation of MILP (compform) formulation
bool useReduction1=true;
bool useReduction2=true;
bool bool_maxC;

// Variables for Reduction
bool Cons_Real_idxs[NB + 1];
int  Real_idxs[NB+1];
bool cons[NB + 1];

// Optimization statistics
int ncall_unique_sort = 0; // Total number of quick sort calls.int ncall_tot = 0;
int ncall = 0;             // Number of quick sort calls per iteration
int nintervals_tot=0;      // Total number of intervals [p^k,p^(k+1)], i.e., total number of KP linear relaxations that has been computed
int nb_ruled_out_tot=0;    // Total number of ruled out items.
double closethegap_time;
double cutoffvalue;        // Best LB found

double blb_before_gap;
double bub_before_gap;
double bounding_proc_time;


// Read input file
int readInput(bool prt) {
    FILE *in;
    int idummy;
    prt=false;
    in = fopen(filename, "r");
    if (in==NULL){
        printf(" *** Warning: cannot open input file %s \n", filename);
        return 1;
    }

    fscanf(in, "%i", &n);
    if (prt) printf("\n Nr. Items %4i \n\n", n);

    if (prt) printf("--i-- -----v----- ---p--- ---c--- \n");
    fscanf(in, "%i %lf", &idummy, &vd[0]);
    vi[0] = int(vd[0] * 1000000 + 0.5);
    if (prt) printf(" %4i %11.8lf       -       -\n", 0, vd[0]);
    for (int i = 1; i <= n; ++i) {
        fscanf(in, "%i %lf %lf %lf", &idummy, &vd[i], &p[i], &c[i]);
        vi[i] = int(vd[i] * 1000000 + 0.5);
        if (prt) {
            printf(" %4i %11.8lf %7.2lf %7.2lf \n", i, vd[i], p[i], c[i]);
            if ((i % 20 == 0) && (i < n))
                printf("--i-- -----v----- ---p--- ---c--- \n");
        }
    }
    fclose(in);

    printf(" Read input file \n");
    return 0;
}

// Generate a set of benchmark instances
int instGenSet(int p_ninst, int p_n, double p_phi, double p_gamma, double p_prec, bool prt) {
    FILE *out;
    char name[100], name2[100], strn[10], sn[20], sphi[20], sgamma[20], sinst[20], txt[5];

    double rx[NB + 1], rxt;
    int vit;

    printf("\n Generating %i instances with %i products [phi %4.2f - gamma %4.2f] \n\n", p_ninst, p_n, p_phi, p_gamma);

    srand(0);

    if (p_n > NB) {
        printf(" *** Warning: cannot handle %i products --> increase NB \n", p_n);
        return 1;
    }
    n = p_n;
    char prefix[50];
    strcpy(prefix, "../MNL_instances/");

    for (int iinst = 1; iinst <= p_ninst; ++iinst) {
        //_itoa_s(p_n, strn, 10);
        sprintf(strn,"%d",p_n);
        sn[0] = 'n';
        sn[1] = '\0';
        strcat(sn, strn);

        sphi[0] = 'p';
        if ((p_phi < 0.76) && (p_phi > 0.74)) {
            sphi[1] = '7';
            sphi[2] = '5';
            sphi[3] = '\0';
        } else if ((p_phi < 0.26) && (p_phi > 0.24)) {
            sphi[1] = '2';
            sphi[2] = '5';
            sphi[3] = '\0';
        } else {
            printf(" *** Warning: cannot handle this value of phi %lf \n", p_phi);
            return 1;
        }

        sgamma[0] = 'g';
        if ((p_gamma < 1.01) && (p_gamma > 0.99)) {
            sgamma[1] = '1';
            sgamma[2] = '0';
            sgamma[3] = '-';
            sgamma[4] = '\0';
        } else if ((p_gamma < 0.51) && (p_gamma > 0.49)) {
            sgamma[1] = '0';
            sgamma[2] = '5';
            sgamma[3] = '-';
            sgamma[4] = '\0';
        } else if ((p_gamma < 2.01) && (p_gamma > 1.99)) {
            sgamma[1] = '2';
            sgamma[2] = '0';
            sgamma[3] = '-';
            sgamma[4] = '\0';
        } else {
            printf(" *** Warning: cannot handle this value of gamma %lf \n", p_gamma);
            return 1;
        }

        //itoa(iinst, sinst, 10);
        sprintf(sinst,"%d",iinst);

        txt[0] = '.';
        txt[1] = 't';
        txt[2] = 'x';
        txt[3] = 't';
        txt[4] = '\0';

        *name = NULL;
        strcat(name, prefix);
        strcat(name, sn);
        strcat(name, sphi);
        strcat(name, sgamma);
        strcat(name, sinst);
        strcat(name, txt);

        if (prt) printf(" Inst.%3i - %s \n\n", iinst, name);

        //if ((err = fopen_s(&out, name, "w")) != 0) {
        if ((out = fopen(name, "w")) == NULL) {
            printf(" *** Warning: cannot open output file %s\n", name);
            return 1;
        }

        //if ((err = fopen_s(&out, name, "w")) != 0) {
        out = fopen(name, "w");
        if (out == NULL) {
            printf(" *** Warning: cannot open output file %s\n", name);
            return 1;
        }

        if (prt) printf(" Random Numbers: \n");
        rxt = 0.0;
        for (int i = 1; i <= p_n; ++i) {
            rx[i] = (double) rand() / (double) RAND_MAX; //see, e.g., https://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1/6219525
            rxt += rx[i];
            if (prt) {
                printf(" %4i %8.6lf\t", i, rx[i]);
                if (i % 5 == 0)
                    printf("\n");
            }
        }
        if (prt)
            printf(" Sum of Random Numbers %8.6lf \n\n", rxt);

        if (prt) printf(" Preference Weights: \n");
        vit = 0;
        for (int i = 1; i <= n; ++i) {
            vi[i] = int(rx[i] / rxt * p_prec);
            vit += vi[i];
        }
        // Make sure that there are enough non-zero preferences
        while (vit < p_prec) {
            int my = int(rand() % n) + 1;
            ++vi[my];
            ++vit;
        }
        for (int i = 1; i <= n; ++i) {
            vd[i] = vi[i] / p_prec;
            assert(vd[i]>0);
            if (prt) {
                printf(" %4i %8i\t", i, vi[i]);
                if (i % 5 == 0)
                    printf("\n");
            }
        }
        if (prt)
            printf(" Sum of Preference Weights %8i \n", vit);
        vi[0] = int(p_phi * vit / (1 - p_phi));
        vd[0] = vi[0] / p_prec;
        if (prt) printf(" No-Purchase Preference    %8i \n\n", vi[0]);

        if (prt) printf(" Product Margin Profits: \n");
        p[0] = 0;
        for (int i = 1; i <= n; ++i) {
            //p[i] = double((rand() * 1000 + rand() % 1000) % int(MAXP + 1));
            p[i] = (int)((double)rand() / ((double)RAND_MAX) * MAXP);
            p[i] = MAX(p[i], 1.0);
            if (prt) {
                printf(" %4i %8.2lf\t", i, p[i]);
                if (i % 5 == 0)
                    printf("\n");
            }
        }

        if (prt) printf("\n Product Fixed Costs: \n");
        c[0] = 0.0;
        for (int i = 1; i <= n; ++i) {
            if (n==500 && p_phi==0.25 && p_gamma==1.0 && iinst==10 && i==256){
                printf("That's it");
            }
            double val = p_gamma * p[i] * vd[i] / (vd[0] + vd[i]); //added for debug purpose
            c[i] = (double) rand() / (double) RAND_MAX * (p_gamma * p[i] * vd[i] / (vd[0] + vd[i]));
            //c[i] = MAX(c[i], 0.01);
            if (c[i]<0 || c[i]>val) {
                printf("Ciaone\n");
                getchar();
            }
            //printf("%f %f %f\n", vd[i], c[i], val);
            assert (vd[i]>0);
            assert (c[i]>0);
            assert (c[i] < val);
            assert (c[i] < val/p_gamma);
            if (prt) {
                printf(" %4i %7.2lf\t", i, c[i]);
                if (i % 5 == 0)
                    printf("\n");
            }

        }
        printf("\n Instance nr.%i Generated \n\n", iinst);

        fprintf(out, "%4i", n);
        fprintf(out, "\n%4i\t%10.8lf", 0, vi[0] / p_prec);
        for (int i = 1; i <= n; ++i)
            fprintf(out, "\n%4i\t%10.8lf\t%7.2lf\t%7.8lf", i, vi[i] / p_prec, p[i], c[i]);
        fclose(out);


        txt[0] = '.';
        txt[1] = 'b';
        txt[2] = 'a';
        txt[3] = 't';
        txt[4] = '\0';

        *name2 = NULL;
        strcat(name2, sn);
        strcat(name2, sphi);
        strcat(name2, sgamma);
        strcat(name2, sinst);
        strcat(name2, txt);

//        if ((err = fopen_s(&out, name2, "w")) != 0) {
        out = fopen(name2, "w");
        if (out == NULL) {
            printf(" *** Warning: cannot open output file \n");
            return 1;
        }

        fprintf(out, "copy %s ..\\\\input.txt", name);
        fclose(out);
    }

    return 0;
}

// Check Cplex output code
int checkStatus(CPXENVptr env, int status) {
    char errmsg[1024];

    if (status == 0) return 0;
    printf("\n CPLEX error \n");
    CPXgeterrorstring(env, status, errmsg);
    printf(" %s \n", errmsg);
    return status;
}

// Add coefficient Cplex constraint
inline void addCoefConstr(int *matind, double *matval, int *nzcnt, int ind, double coef) {
    if (ind == IMAXI) return;
    matind[*nzcnt] = ind;
    matval[*nzcnt] = coef;
    *nzcnt = *nzcnt + 1;
}

// Compute the profit of a solution/assortment passed through x[] (Cplex vars)
double clcProfitSol(double *x, int *ind, bool prt) {
    // called  in CompForm, lazycallbacks, closeNoGood.
    // x: vector of Cplex Variables
    // ind: ind[j] is the index of item j in the cplex model. x[ind[j]]=1 if j in the assortment, 0 othw.
    // prt: boolean, to print statistics to output.
    double ct, pt, vdt;
    int cnt;
    ct = pt = 0.0;
    cnt = 0;
    vdt = vd[0];
    if (prt)
        printf(" Products: \n");
    for (int j = 1; j <= n; ++j)
        // if j is in the assortment (eq. to x[ind[j]]==1)
        if (x[ind[j]] > 0.9) {
            ct += c[j];
            vdt += vd[j];
            ++cnt;
        }

    for (int j = 1; j <= n; ++j)
        if (x[ind[j]] > 0.9)
            pt += p[j] * vd[j];
    pt /= vdt;

    if (prt) {
        printf("\n Sum of preferences        %10.8lf \n", vdt);
        printf(" Total costs           %10.4lf \n", ct);
        printf(" Total marginal profit %10.4lf \n", pt);
        printf(" Final profit          %10.8lf \n", pt - ct);
        printf(" Items in solution     %4i \n", cnt);
    }
    return pt - ct;
}

// Compute the profit of a solution/assortment containing the first indicik items in *indici.
double clcProfitSol2(int *indici, int indicik) {
    // Called only in parametricUB
    double ct, pt, vdt;
    int jj, j;

    ct = pt = 0.0;
    vdt = vd[0];
    for (jj = 1; jj <= indicik; ++jj) {
        j = indici[jj];
        ct += c[j];
        vdt += vd[j];
        pt += p[j] * vd[j];
    }
    pt /= vdt;
    return pt - ct;
}

// (MILP) Compact formulation Kunnumkal and Mart�nez-de-Alb�niz, Oper Res (2019)
double compForm(bool *inside, double tlim, int mcard) {
    /*
     * Single class
     */
    double st = clock();
    CPXENVptr env;
    CPXLPptr lp;
    int status;
    double obj[1], rhs[1], lb[1], ub[1], matval[2 * NB + 2];
    int matind[2 * NB + 2], matbeg[1], nzcnt;
    char sense[1], xctype[1], name[50], *colname[1], *rowname[1];
    double objval, x[2 * NB + 3];
    int ind_u[NB + 1], ind_x[NB + 1], nvar;

    for (int i=0; i<=n; ++i) Real_idxs[i]=i;

    printf("\n\n Compact Formulation \n\n");

    objval = 0.0;

    env = CPXopenCPLEX(&status);
    if (checkStatus(env, status)) return objval;

    lp = CPXcreateprob(env, &status, "AO");
    if (checkStatus(env, status)) return objval;

    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
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

    colname[0] = name;

    nvar = 0;

    //  variables u, and their obj coefficients
    lb[0] = 0.0;
    xctype[0] = 'C';
    for (int j = 0; j <= n; ++j) {
        obj[0] = p[j];

        sprintf_s(name, sizeof(char) * 12, "u%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, NULL, xctype, colname);
        if (checkStatus(env, status)) return objval;
        ind_u[j] = nvar++;
    }

    // variables x, and their obj coefficients
    xctype[0] = 'B';
    lb[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        obj[0] = -c[j];
        ub[0] = 1.0;
        if (!inside[j]) ub[0] = 0.0;

        sprintf_s(name, sizeof(char) * 12, "x%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
        if (checkStatus(env, status)) return objval;
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
        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // Ctr (2.c) in the paper
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], vd[0] + vd[j]);
        addCoefConstr(matind, matval, &nzcnt, ind_x[j], -vd[j]);
        sprintf_s(name, sizeof(char) * 15, "Rel.u%dx%d", j, j);
        status = CPXaddrows(env, lp, 0, 1, 2, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // Ctr (2.d) in the paper
    sense[0] = 'E';
    rhs[0] = 1.0;
    nzcnt = 0;
    for (int j = 0; j <= n; ++j)
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], 1.0);
    sprintf_s(name, sizeof(char) * 15, "Sum.u");
    status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
    if (checkStatus(env, status)) return objval;

    // Cardinality constraint
    if (mcard < n) {
        sense[0] = 'L';
        rhs[0] = mcard;
        nzcnt = 0;
        for (int j = 1; j <= n; ++j)
            addCoefConstr(matind, matval, &nzcnt, ind_x[j], 1.0);
        sprintf_s(name, sizeof(char) * 15, "Card");
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    if (checkStatus(env, status)) return objval;

    // Solve Linear Relaxation
    if (solveLP){
        CPXchgprobtype(env, lp, 0);
        status = CPXlpopt(env,lp);
        // Opt Val
        status = CPXgetobjval(env, lp, &objval);
        double et = clock();
        fprintf(csv_out, "%f, %f\n", objval, (et-st)/(double)CLOCKS_PER_SEC);
        printf("Objective Value: %f", objval);
        return objval;
    }

    // Solve MIP
    status = CPXmipopt(env, lp);

    if (checkStatus(env, status)) return objval;

    // Get Opt Sol Stats:
    // Termination code (optimal vs. time limit vs optimal within tolerance)
    int ex_code= CPXgetstat(env, lp);

    double gap;
    CPXgetmiprelgap(env, lp, &gap);

    // Opt Obj Val
    status = CPXgetobjval(env, lp, &objval);
    if (checkStatus(env, status)) return objval;
    double et = clock(); //mine

    double dual_bound;
    status = CPXgetbestobjval(env, lp, &dual_bound);

    int N_nodes =  CPXgetnodecnt(env, lp);

    printf("\n Optimal solution value %12.4f \n\n", objval);
    // Save stats to file
    fprintf(csv_out, "%f, %d, %f, %f, %f,%d\n", objval, ex_code, dual_bound, gap,(et-st)/(double)CLOCKS_PER_SEC, N_nodes);

    // Get Optimal solution
    status = CPXgetx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1);
    if (checkStatus(env, status)) return objval;

    // Manually compute profit, and print solution
    clcProfitSol(x, ind_x, true);
    printf("Optimal Solution from CompForm:\n");
    for (int j = 1; j <= n; ++j)
        if (x[ind_x[j]] > 0.9) {
            printf("%d\n", j);
        }
    status = CPXfreeprob(env, &lp);
    if (checkStatus(env, status)) return objval;

    status = CPXcloseCPLEX(&env);
    if (checkStatus(env, status)) return objval;

    return objval;
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

    if (l < r) {
        j = QS_Partition(vett, indici, l, r);
        QS_QuickSort(vett, indici, l, j - 1);
        QS_QuickSort(vett, indici, j + 1, r);
    }
}

// Solve KP Linear Relaxation for every interval
double LPRelKPbz(double *pr, double *w, double cap, int *indici, double oldrapp, bool *taken, double *minr) {
    double tw, ub, r[NB + 1], maxr;
    int maxp, j;

    tw = ub = 0.0;
    *minr = IMAXI;
    for (j = 1; j <= n; ++j) {
        taken[j] = false;
        r[j] = pr[j] / w[j];
        // Insert in the solution all elements j: r[j]>= oldrapp
        if (r[j] >= oldrapp) {
            tw += w[j];
            ub += pr[j];
            *minr = pr[j] / w[j];//cs: MIN(*minr, r[j])
            taken[j] = true;
        }
    }

    // if too many items (e.g., first iteration) --> solve the KP relaxation from scratch.
    if (tw > cap) {
        // remove ALL elements from the KP
        for (j = 1; j <= n; ++j)
            taken[j] = false;

        // sort by p/w
        ++ncall_unique_sort; //stats counting how many times we need to perform a quicksort
        ++ncall;
        QS_QuickSort(r, indici, 1, n);

        ub = tw = 0.0;

        // Insert until critical item
        for (int j = 1; j <= n; ++j) {
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
        // Fill remainining capacity
        int my_cnt = 0;
        while (tw < cap) {
            ++my_cnt;
            maxp = 0;
            maxr = 0.0;
            for (j = 1; j <= n; ++j)
                if (!taken[j])
                    if (maxr < r[j]) {
                        maxr = r[j];
                        maxp = j;
                    }
            if (maxp == 0) {
                break;
            }

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
    }

    return ub;
}

// Solve the KP Linear Relaxation with Cardinality constraint
double LPRelKPbzcard(double *pr, double *w, double cap, int *indici, double oldrapp, int olditer, bool *taken, double *minr,
              int *lastiter) {
    // oldrapp: ratio of the critical item p_s/w_s at the previous iteration.
    // olditer: to speedup convergence, use the number of iters performed to update lambda in the previous iteration,
    //          and start from this set of updates with lambda = lambdab * olditer / double(iterb);
    // lastiter: pointer to save the number of iter/lambdab updates used in this iteration. This will be used as starting
    //            olditer in the next iteration.
    double tw, ub, r[NB + 1], bub, lambdab, maxr;
    int maxp, j, counter, iterb, iter;
    bool takenb[NB + 1]; // Solution associated to best UB, while take is the sol associated to each iter.

    // Initial guess for lambda
    lambdab = 0.0;
    for (int j = 1; j <= n; ++j)
        lambdab = MAX(lambdab, pr[j]);
    iterb = 100000;

    bub = IMAXI;
    // Improve lambda iteratively
    for (iter = MAX(0, olditer); iter <= iterb; ++iter) {
        double lambda = lambdab * iter / double(iterb);
        double prl[NB + 1];
        int indicil[NB + 1];

        // Compute ratios
        for (int j = 1; j <= n; ++j) {
            prl[j] = pr[j] - lambda; // Add lambda to the item-cost
            r[j] = prl[j] / w[j];
            indicil[j] = indici[j];
        }

        for (j = 1; j <= n; ++j)
            taken[j] = false;

        ub = lambda * maxC; // Obj of Cardinality constrained Lagrangian Relax
        tw = 0.0;
        counter = 0;
        // Fill with items based on oldrapp
        for (int j = 1; j <= n; ++j)
            if (r[j] >= oldrapp) {
                ub += prl[j];
                tw += w[j];
                taken[j] = true;
                *minr = MIN(*minr, prl[j] / w[j]);
                ++counter;
            }

        if (tw <= cap) {
            // Fill remaining capacity
            while (tw < cap) {
                maxp = 0;
                maxr = 0.0;
                for (int j = 1; j <= n; ++j)
                    if (!taken[j])
                        if (maxr < r[j]) {
                            maxr = r[j];
                            maxp = j;
                        }
                if (maxp == 0) {
                    break;
                }

                if (w[maxp] <= cap - tw) {
                    ub += prl[maxp];
                    tw += w[maxp];
                    taken[maxp] = true;
                    *minr = prl[maxp] / w[maxp];
                    ++counter;
                    if (cap - tw < 0.00000000001)
                        break;
                } else {
                    ub += (cap - tw) / w[maxp] * prl[maxp];
                    tw = cap;
                    *minr = prl[maxp] / w[maxp];
                    ++counter;
                    break;
                }
            }
            goto CHECKUB;
        }

        // Sort item and do the relax from scratch
        ++ncall_unique_sort;
        ++ncall;
        QS_QuickSort(r, indicil, 1, n);

        for (j = 1; j <= n; ++j)
            taken[j] = false;

        ub = lambda * maxC; // Obj of Cardinality constrained Lagrangian Relax
        tw = 0.0;
        counter = 0;
        for (int j = 1; j <= n; ++j) {
            maxp = indicil[j];
            if (prl[maxp] / w[maxp] <= 0.0)
                break;

            if (w[maxp] <= cap - tw) {
                ub += prl[maxp];
                tw += w[maxp];
                taken[maxp] = true;
                *minr = prl[maxp] / w[maxp];
                ++counter;
                if (cap - tw < 0.00000000001) break;
            } else {
                // Insert critical item
                ub += (cap - tw) / w[maxp] * prl[maxp];
                tw = cap;
                *minr = prl[maxp] / w[maxp];
                ++counter;
                break;
            }
        }

        CHECKUB:;
        oldrapp = MIN(oldrapp, *minr);
        if (ub < bub) {
            bub = ub;
            for (int j = 1; j <= n; ++j)
                takenb[j] = taken[j];
        } else {
            break;
        }

        if (counter <= maxC) {
            break;
        }
    }
    if (iter>iterb) iter=iterb;
    *lastiter = iter;

    for (int j = 1; j <= n; ++j)
        taken[j] = takenb[j];

    return bub;
}

// (MILP+) 
double closethegap_compform(double tlim, double pL, double pU, double p_lb) {
    /*
     * Same formulation of CompForm, but with
     */
    double st = clock();
    CPXENVptr env;
    CPXLPptr lp;
    int status;
    double obj[1], rhs[1], lb[1], ub[1], matval[2 * NB + 2];
    int matind[2 * NB + 2], matbeg[1], nzcnt;
    char sense[1], xctype[1], name[50], *colname[1], *rowname[1];
    double objval, x[2 * NB + 3];
    int ind_u[NB + 1], ind_x[NB + 1], nvar;

    printf("\n\n Compact Formulation \n\n");

    objval = 0.0;

    env = CPXopenCPLEX(&status);
    if (checkStatus(env, status)) return objval;

    lp = CPXcreateprob(env, &status, "AO");
    if (checkStatus(env, status)) return objval;

    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
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
        if (checkStatus(env, status)) return objval;
        ind_u[j] = nvar++;
    }

    // variables x, and their obj coefficients
    xctype[0] = 'B';
    lb[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        obj[0] = -c[j];
        if (Cons_Real_idxs[Real_idxs[j]]){
            ub[0] = 1.0;
        }
        else {
            ub[0] = 0.0; // Ruled out items
            assert (!Cons_Real_idxs[Real_idxs[j]]);
        }

        sprintf_s(name, sizeof(char) * 12, "x%d", j);
        status = CPXnewcols(env, lp, 1, obj, lb, ub, xctype, colname);
        if (checkStatus(env, status)) return objval;
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
        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // Ctr (2.c) in the paper
    sense[0] = 'L';
    rhs[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], vd[0] + vd[j]);
        addCoefConstr(matind, matval, &nzcnt, ind_x[j], -vd[j]);
        sprintf_s(name, sizeof(char) * 15, "Rel.u%dx%d", j, j);
        status = CPXaddrows(env, lp, 0, 1, 2, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
    }

    // Ctr (2.d) in the paper
    sense[0] = 'E';
    rhs[0] = 1.0;
    nzcnt = 0;
    for (int j = 0; j <= n; ++j)
        addCoefConstr(matind, matval, &nzcnt, ind_u[j], 1.0);
    sprintf_s(name, sizeof(char) * 15, "Sum.u");
    status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
    if (checkStatus(env, status)) return objval;

    // MILP+, ctr (12) in the paper, <=
    sense[0] = 'L';
    rhs[0] = pU;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], 1/vd[0]);
        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
        break;
    }

    // MILP+, ctr (12) in the paper, >=
    sense[0] = 'G';
    rhs[0] = pL;
    for (int j = 1; j <= n; ++j) {
        nzcnt = 0;
        addCoefConstr(matind, matval, &nzcnt, ind_u[0], 1/vd[0]);
        sprintf_s(name, sizeof(char) * 15, "Rel.u%du0", j);
        status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rowname);
        if (checkStatus(env, status)) return objval;
        break;
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
        if (checkStatus(env, status)) return objval;
    }

    if (checkStatus(env, status)) return objval;

    // Solve Linear Relaxation
    if (solveLP){
        CPXchgprobtype(env, lp, 0);
        status = CPXlpopt(env,lp);
        // Opt Val
        status = CPXgetobjval(env, lp, &objval);
        double et = clock();
        fprintf(csv_out, "%f, %f\n", objval, (et-st)/(double)CLOCKS_PER_SEC);
        printf("Objective Value: %f", objval);
        return objval;
    }

    // Solve MIP
    status = CPXmipopt(env, lp);

    if (checkStatus(env, status)) return objval;

    // Get Opt Sol Stats:
    // Termination code (optimal vs. time limit vs optimal within tolerance)
    int ex_code= CPXgetstat(env, lp);

    double gap;
    CPXgetmiprelgap(env, lp, &gap);

    // Opt Obj Val
    status = CPXgetobjval(env, lp, &objval);
    if (checkStatus(env, status)) return objval;
    double et = clock(); //mine

    double dual_bound;
    status = CPXgetbestobjval(env, lp, &dual_bound);

    int N_nodes =  CPXgetnodecnt(env, lp);

    printf("\n Optimal solution value %12.8f \n\n", objval);

    // Get Optimal solution
    status = CPXgetx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1);
    if (checkStatus(env, status)) return objval;

    // Manually compute profit, and print solution
    double manual_prof = clcProfitSol(x, ind_x, true);

    status = CPXfreeprob(env, &lp);
    if (checkStatus(env, status)) return objval;

    status = CPXcloseCPLEX(&env);
    if (checkStatus(env, status)) return objval;

    return manual_prof;
}

double get_K(double t, double rho){
    // t = either tmin or tmax
    double K = log(t) / log(1+rho);
    return K;
}

void parametricUB_topal(double fstep, double lstep) {
    /*
     * fstep: first level of grid density (granularity rho in the paper: how many digits)
     * lstep: last  level of grid density
     * incr = 50: nb of intervals used to split [pmin, pmax] after the last iteration, for the "enumerative phase"
     */
    double *lbs, *ubs, bestlb, bestub, *rappmin;
    int denomin, denomout, indici[NB+1], maxiterlast;
    double pmax, pmin;
    double profitti[NB + 1], pesi[NB + 1];
    bool taken[NB + 1], wcard; //wcard: wether cardinality ctr is not satisfied.

    // To remember the original item idx, after reductions (and compare solutions with other methods)
    for (int i=0; i<=n; ++i) {
        Real_idxs[i]=i;
        if (i!=0) Cons_Real_idxs[i] = true;
    }

    double st, et, bounding_start;
    double rho = fstep;
    st = clock();
    bounding_start = clock();
    printf("\n\n Parametric Upper Bound \n\n");

    bestlb = 0.0;
    wcard = true;

    // Compute pmin, pmax
    double pref_sum = 0;
    double pref_min = IMAXI;
    for (int i =1; i <= n; ++i) {
        pref_sum += vd[i];
        if (vd[i]>0) pref_min = MIN(pref_min, vd[i]);
    }

    // Starting range, as in Feldman and Topaloglu (2015) and Kunnumkal (2019) (Unnormalied v0, as in the latter)
    pmin = 1/(vd[0] + pref_sum);
    pmax = 1/(vd[0] + pref_min);

    AGAIN:; // for every rho in [fstep, lstep]:

    // Get exponents K
    int k_min = int(get_K(pmin, rho));
    if (k_min<0) k_min--; //Should always be the case.
    int k_max = int(get_K(pmax, rho));
    if (k_max>0) k_max++;
    assert (pow(1+rho, k_min) <= pmin);
    assert (pow(1+rho, k_max) >= pmax - 1e-8);

    printf(" Rho = %10.8lf \n", rho);
    maxiterlast = 0;


    printf("k_min: %d  k_max:%d\n", k_min, k_max);

    // Keep track over lower/upper bound for each denomitor (useful for "pruning" the [pmin,pmax] range)
    lbs = (double *) malloc((k_max - k_min + 1) * sizeof(double));
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
        if ((maxC >= n) || (!wcard)) {
            double newub;
            newub = LPRelKPbz(profitti, pesi, 1/pL - vd[0], indici, oldrapp, taken, &minr);

            ubs[k-1-k_min] = newub;

            rappmin[k-1-k_min] = minr;
        }
        else {
            double ubb;

            if (k == k_max) olditer = 0;
            else olditer = newiter;

            ubb = LPRelKPbzcard(profitti, pesi, 1/pL-vd[0], indici, oldrapp, olditer, taken, &minr, &newiter);
            ubs[k-1-k_min] = ubb;
            maxiterlast = MAX(maxiterlast, newiter);
            rappmin[k-1-k_min] = minr;
        }

        int mycard = 0;
        double tp = 0;
        int insol[NB + 1], insolk;

        // compute insol (items in the assortment)
        if (maxC >= n) {
            insolk = 0;
            // obtain LowerBound from the solution of the linear relaxation (removing critical item)
            for (int j = 1; j <= n; ++j) {
                if (taken[j]) {
                    tp += vd[j];

                    ++insolk;
                    insol[insolk] = j;
                }
            }
        } else {
            insolk = 0;

            bool taken2[NB + 1];

            for (int j = 1; j <= n; ++j)
                taken2[j] = false;

            mycard = 0;
            tp = 0;

            // As for the unconstrained case, get solution from LP relaxation
            for (int j = 1; j <= n; ++j)
                if (taken[j]) {
                    tp += vd[j];
                    ++mycard;
                }

            // If the cardinality ctr. is violated
            while (mycard > maxC) {
                int    jmin = 0;
                double vmin = IMAXI;

                // Find item of minimum profit
                for (int j = 1; j <= n; ++j)
                    if (taken[j] && !taken2[j])
                        if (vd[j] < vmin) {
                            vmin = vd[j];
                            jmin = j;
                        }
                if (jmin == 0) {
                    insolk = 0;
                    break;
                }
                // Remove it from the knapsack
                tp -= vd[jmin];
                --mycard;
                taken2[jmin] = true;
            }

            // Build feasible solution (for lb computation)
            for (int j = 1; j <= n; ++j)
                if (taken[j] && !taken2[j]) {
                    ++insolk;
                    insol[insolk] = j;
                }
        }

        // Update bounds
        if (ubs[k-1-k_min] > bestlb) {
            lbs[k-1-k_min] = clcProfitSol2(insol, insolk);
        } else
            lbs[k-1-k_min] = 0.0;

        bestub = MAX(bestub, ubs[k-1-k_min]);
        bestlb = MAX(bestlb, lbs[k-1-k_min]);
    }
    printf(" BestLB %lf - BestUB %lf \n", bestlb, bestub);
    et = clock();
    printf(" CPU time Param UB %lf \n", (et - st) * 0.001);
    printf(" #calls QuickSort %i - [%i] \n", ncall, k_max - k_min);
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
                pmin = pL; // cs: minin= smallest t: z(t, pho)> bestlb
                new_kmin= k;
            }
            pmax = pU; // cs: maxin= biggest t: z(t, pho)> bestlb
            new_kmax= k+1;//k+1
        }
    }
    k_min = new_kmin;
    k_max = new_kmax;

    printf(" Intervals %i - In %i - Out %i - Range %lf-%lf \n", k_max - k_min, denomin, denomout, pmin, pmax);

    for (int jj = 0; jj <= n; ++jj)
        cons[jj] = true;

    // Counter of ruled out items
    int cntout = 0;

    // Reduction procedures
    if ((!useReduction1)) goto SKIPRED1;

    // First reduction
    st = clock();
    printf("\n Ruling out items: \n");
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
    printf("\n Ruled out %i [%i] \n", cntout, n);
    et = clock();
    printf(" CPU time Ruling Out Items %lf \n", (et - st) * 0.001);

    SKIPRED1:;
    if ((!useReduction2)) goto SKIPRED;

    // Second Reduction procedure
    st = clock();
    printf("\n Ruling out items: ");
    cntout = 0;
    for (int j = 1; j <= n; ++j) {
        if (!cons[j]) continue;

        for (int k = k_min; k <= k_max-1; ++k) {
            double prof;
            double pL= pow(1+rho, k);
            double pU= pow(1+rho, k+1);
            prof =  pU * p[j] * vd[j] - c[j];

            // Compute the bound only for products such that x_j=0 when computing the linear relax z(t, rho)
            if (prof / vd[j] >= rappmin[(k)-k_min])
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
    printf("\n Ruled out %i [%i] \n", cntout, n);
    et = clock();
    printf(" CPU time Ruling Out Items %lf \n", (et - st) * 0.001);

    SKIPRED:;
    int counter = 0;
    // Update problem variables
    for (int j = 1; j <= n; ++j) {
        if (cons[j]) {
            ++counter;
            p[counter] = p[j];
            c[counter] = c[j];
            vd[counter] = vd[j];
            vi[counter] = vi[j];
            Real_idxs[counter] = Real_idxs[j];
        }
    }
    printf(" Problem reduced to %i items (from %i) \n\n", counter, n);
    n = counter;

    if (rho > lstep+1e-8) {
        if (maxiterlast == 0)
            wcard = false;
        else
            wcard = true;

        free(lbs);
        free(ubs);
        free(rappmin);

        rho /= 10.0;
        goto AGAIN;
    }
    printf("\n\n Closing the Gap \n\n");
    cutoffvalue = bestlb;
    blb_before_gap = bestlb;
    bub_before_gap = bestub;
    double bounding_end = clock();
    bounding_proc_time = (bounding_end - bounding_start)/CLOCKS_PER_SEC;

    // Solve MILP+ for closing the gap
    st = clock();
    cutoffvalue = closethegap_compform(600., pmin, pmax, cutoffvalue);
    closethegap_time = (clock()-st)/(double)CLOCKS_PER_SEC;


    printf("Optimal cutoff:                %lf\n", cutoffvalue);
    free(lbs);
    free(ubs);
    free(rappmin);

    printf("Solution with Original Variable indexes:\n");
    printf("Old mapping:\n");
}

int main(int argc, char* argv[]) {
    double st, et;
    int inst_type = 0;
    bool inside[NB + 1];
    double lstep=1e-7;

    /*
    // Uncomment this to generate instances. Parameters:
    // 1. Number of instances
    // 2. Number of items
    // 3. Parameter (phi in the paper) determining the probability of the no-purchase option.
    // 4. Parameter (Gamma in the paper) determining how big fixed costs are. (bigger gamma -> bigger product costs)
    // 5. Preicion: number of digits of the product preferences.
    // 6. boolean parameter to print some data generation statistics.
    */
//    instGenSet(50, 100, 0.25, 0.5, 100000000.0, false);
//    instGenSet(50, 100, 0.25, 1.0, 100000000.0, false);
//    instGenSet(50, 100, 0.75, 0.5, 100000000.0, false);
//    instGenSet(50, 100, 0.75, 1.0, 100000000.0, false);
//    instGenSet(50, 200, 0.25, 0.5, 100000000.0, false);
//    instGenSet(50, 200, 0.25, 1.0, 100000000.0, false);
//    instGenSet(50, 200, 0.75, 0.5, 100000000.0, false);
//    instGenSet(50, 200, 0.75, 1.0, 100000000.0, false);
//    instGenSet(50, 500, 0.25, 0.5, 100000000.0, false);
//    instGenSet(50, 500, 0.25, 1.0, 100000000.0, false);
//    instGenSet(50, 500, 0.75, 0.5, 100000000.0, false);
//    instGenSet(50, 500, 0.75, 1.0, 100000000.0, false);
//    instGenSet(50, 1000, 0.25, 0.5, 100000000.0, false);
//    instGenSet(50, 1000, 0.25, 1.0, 100000000.0, false);
//    instGenSet(50, 1000, 0.75, 0.5, 100000000.0, false);
//    instGenSet(50, 1000, 0.75, 1.0, 100000000.0, false);
//    return 0;
    /* Uncomment above for generating instances */


    // Read Input arguments
    if (argc < 4){
        printf("Error: need the following arguments \n 1.input filename\n 2.maxC={true, false}\n");
        printf(" 3.(bool)mip \n 4.(Optional) id_string for the experiment\n ");
        return -1;
    }

    // full path of input instance
    strcpy(filename, argv[1]);

    // cardinality constrained or not
    bool_maxC = (strcmp(argv[2],"true")==0 | strcmp(argv[2],"True")==0) ? true : false;

    // if 1 --> MILP, if 0--> LLRS
    bool mip = (strcmp(argv[3],"true")==0 | strcmp(argv[3],"True")==0) ? 1 : 0;

    // id_string
    if (argc>=5) strcpy(id_string, argv[4]);
    else         strcpy(id_string, "no_id_string");

    // Output files
    char csv_name[1000];
    if (solveLP) {
        assert (mip);
        strcpy(csv_name, "../output_lp.csv");
    } else if(mip)  strcpy(csv_name,"../output_mip.csv");
    else            strcpy(csv_name,"../output_LLRS.csv");

    if ((csv_out = fopen(csv_name, "a")) == NULL) {
        printf(" *** Warning: cannot open csv file %s\n", csv_name);
        return 1;
    }

    // Read Instances
    if (readInput(false)) return -1;
    fprintf(csv_out, "%s, %s, %d, %d,", id_string, filename, bool_maxC, mip);

    for (int i = 0; i <= n; ++i)
        inside[i] = true;

    // Set cardinality constraint
    if (bool_maxC) maxC = int(n / 2);
    else maxC = n;

    if (mip == 1) {
        // Solve MILP formulation
        double mipub;

        st = clock();
        mipub = compForm(inside, 600, maxC);
        et = clock();

        FILE *out;
        if ((out = fopen("out.txt", "a")) == NULL) {
            printf(" > Cannot open output file \n");
            exit(0);
        } else {
            fprintf(out, "%10.4lf\t%7.2lf\n", mipub, (et - st) * 0.001);
            fclose(out);
        }
        fclose(csv_out);
        return 0;
    }

    // Solve using LLRS
    st = clock();

    parametricUB_topal(1e-2, lstep); // 1.0e7, 1000

    et = clock();

    // Print Statistics (both to screen and to and to file)
    printf(" Best LB %10.6lf \n", cutoffvalue);
    printf(" CPU time %lf \n", (et - st) * 0.001);


    fprintf(csv_out, "%f, %f, %d,", cutoffvalue, (et-st)/(double)CLOCKS_PER_SEC, nintervals_tot);
    fprintf(csv_out, "%d, %d,", useReduction1, useReduction2);
    fprintf(csv_out, "%d,", nb_ruled_out_tot);
    fprintf(csv_out, "%lf, %lf, %lf, %lf\n", closethegap_time, blb_before_gap, bub_before_gap, bounding_proc_time);


    FILE *out;
    if ((out = fopen("out.txt", "a")) == NULL) {
        printf(" > Cannot open output file \n");
        exit(0);
    } else {
        fprintf(out, "%10.4lf\t%7.2lf\n", cutoffvalue, (et - st) * 0.001);
        fclose(out);
    }
    printf("Nb sorting total: %d\n", ncall_unique_sort);
    printf("Nb ruled out tot: %d\n", nb_ruled_out_tot);
    printf("Nb intervals tot: %d\n", nintervals_tot);
    printf("Close the gap time: %lf\n", closethegap_time);
    fclose(csv_out);

    return 0;
}