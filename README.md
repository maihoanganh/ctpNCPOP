# ctpNCPOP
- ctpNCPOP is a Julia package to solve noncommutative polynomial optimization problem (NCPOP), in particular eigenvalue minimization over tuples of symmetric matrices:

**𝜆_min = inf_{A in S^n} { v*f(A)v : gi(A) >= 0, hj(A) = 0, |v|=1 },**

with some special cases of inequality constraints **gi**:

Case 1: Annulus constraints on subsets of variables: **Ui >= ||x(Ti)||^2 >= Li**.

Case 2: Simplex constraints: **xi >= 0, 1 - x1 -...- xn >= 0**.

- The main ingredient of ctpNCPOP is to solve Moment-SOS (semidefinite) relaxations of the form:

**v = inf_X { <C,X> : X is positive semidefinite, AX = b },**

such that the constant trace property (CTP) holds:

**AX = b => trace(X) = a,**

by using Conditional gradient-based augmented Lagrangian (CGAL).

- Although possibly slower than the classical interior-point method on sparse problems, ctpNCPOP is more efficent on the dense ones.

- ctpNCPOP can combine CTP with correlative sparsity (CS), term sparsity (TS) and their combination (CS-TS) to overcome memory issues encountered when considering semidefinite relaxations of large-scale problems with sparse data.


# Required softwares
ctpNCPOP has been tested with the following software:
- Ubuntu 18.04.4
- Julia 1.3.1

The following interior-point semidefinite solver is used for comparison purposes:
- [Mosek 9.1](https://www.mosek.com)

Before installing ctpNCPOP, you should install [TSSOS](https://github.com/wangjie212/NCTSSOS) with the following commands:
```ruby
Pkg> add https://github.com/wangjie212/NCTSSOS
```

# Installation
- To use ctpNCPOP in Julia, run
```ruby
Pkg> add https://github.com/maihoanganh/ctpNCPOP.git
```

# Usage
The following examples briefly guide you to use ctpNCPOP:

## Noncommutative polynomial optimization
Consider the following NCPOP on the unit ball:
```ruby
using DynamicPolynomials

@ncpolyvar x[1:2] # variables

f=x[1]^2+0.5*x[1]*x[2]+x[2]^2 # the objective polynomial to minimize

g=[1.0-sum(x.^2)] # the inequality constraints
h=[x[1]+x[2]-0.5] # the equality constraints

k=2 # relaxation order

using ctpNCPOP

# get information from the input data f,gi,hj
n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f,dg,dh=ctpNCPOP.get_info(x,f,g,h)

# get the optimal value of the Moment-SOS relaxation of order k
opt_val=ctpNCPOP.POP_dense_CGAL(  n, # the number of variables
                                m, # the number of the inequality constraints
                                l, # the number of the equality constraints
                                lmon_g, # the number of terms in each inequality constraint
                                supp_g, # the support of each inequality constraint
                                coe_g, # the coefficients of each inequality constraint
                                lmon_h, # the number of terms in each equality constraint
                                supp_h, # the support of each equality constraint
                                coe_h, # the coefficients of each equality constraint
                                lmon_f, # the number of terms in the objective polynomial
                                supp_f, # the support of the objective polynomial
                                coe_f, # the coefficients of the objective polynomial
                                dg, # the degree of each inequality constraint
                                dh, # the degree of each equality constraint
                                k,
                                maxit=Int64(1e6), # the maximal iteration of CGAL solver
                                tol=1e-3, # the tolerance of CGAL solver
                                use_eqcons_to_get_constant_trace=false, # use the equality constraints to get constant trace
                                check_tol_each_iter=true ) # check the tolerance at each iteration
```

See other examples in this [link](https://github.com/maihoanganh/ctpNCPOP/tree/main/examples).


# References
For more details, please refer to:

**N. H. A. Mai, A. Bhardwaj, V. Magron. The Constant Trace Property in Noncommutative Optimization, 2021. Forthcoming.**

The following codes allow one to run the paper's benchmarks:
```ruby
using ctpNCPOP


ctpNCPOP.test_dense_POP_ball(10,30,2,have_eqcons=true) # Table 1

ctpNCPOP.test_dense_POP_box(10,30,2,have_eqcons=true) # Table 2

ctpNCPOP.test_CS_POP_ball(1000,10,20,2,have_eqcons=true) # Table 3

```
