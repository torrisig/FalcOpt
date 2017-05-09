<img src="https://raw.githubusercontent.com/torrisig/FalcOpt/master/logo/logo.jpg"  width="250"  />

# FalcOpt
First-order Algorithm via Linearization of Constraints for OPTimization

# Solving efficient nonlinear Model Predictive Control problems with FalcOpt

**FalcOpt stands for First-order Algorithm via Linearization of Constraints for OPTimization. It allows the user to solve nonlinear Model Predictive Control problems via a novel first-order algorithm presented in [[1](https://arxiv.org/abs/1610.06834)]. This tool works in Matlab and it allows one to automatically generate the C code to solve a custom-made control problem, together with the MEX interface. It heavily exploits sparsity of the problem and it is specifically indicated for time-critical applications.**

The nonlinear Model Predictive Control problem that we aim to solve is:

![<opt problem>](http://mathurl.com/lmdwk5e.png)

where the index _j_ spans the predicted horizon state _x_ and input _u_ in the horizon _N_. The bounds satisfy _a_ < _b_  componentwise, _ck_ > 0 and _Q, P, R_ are positive semidefinite. The discrete-time dynamics, _f_, can be nonlinear and the nonlinear function _n_ can be nonconvex. Besides the dynamics, all of the other constraints can be missing.

## Installation

Add the folder containing +FalcOpt to your Matlab path without subfolders. 

## Setting the problem

To solve the optimization problem presented above, first we run an offline code generation procedure which builds the required C (and possibly mex) functions. Then, the code is called in real-time at every time step _k_ in a receding horizon fashion.

### Code generation

The code is generated via the following command:

`info = falcopt.generateCode(dynamics, N, nx, nu, Q, P, R);`

where:
* `dynamics` is a Matlab function handle containing the nonlinear discretized dynamics _f_;
* `N` is the prediction horizon
* `nx` is the state dimension
* `nu` is the input dimension
* `Q` is the state cost weight
* `R` is the input cost weight
* `P` is the terminal cost weight
* `info` is a Matlab structure containing some general information about the problem (the expected number of FLOPS and the name or header of possible additional files).

The command `falcopt.generateCode` can be executed by specifying additional arguments, which allow one to consider e.g. additional constraints and their automatic differentiation, tolerances or other preferences. The syntax is:

`info = falcopt.generateCode(dynamics, N, nx, nu, Q, P, R, <option>, <value>, ..., <option>, <value>);` 

where the additional options are given in Matlab by executing the command `help falcopt.generateCode`.

### Call the solver

Once that the solver is generated, the function is called in Matlab via the following command:

`[u, flag, info] = my_code(x0, u_sequence);`

where:
* 'my_code' or an user defined name is the mex function calling the homonym C file  
* `x0` is the initial state
* `u_sequence` is a warm-start sequence of the _N_ predicted inputs
and additional arguments such as desired real-time references, terminal constraints parameters can be set in the code generation phase (check out Matlab `help falcopt.generateCode`).

The outputs are the following:
* `u` is the first input, result of the optimization problem, ready to be applied to the system
* `flag` is the exitflag of the optimization problem (>0 successful, 0 max number of iteration reached, <0 an error has occurred - the problem may be infeasible)
* `info` is a struct containing additional information, such as iterations, solve time, optimal value of the objective function and full predicted states and inputs. The amount of information can be selected during the code generation phase (check out Matlab `help falcopt.generateCode`). 

## FAQ

**What compiler is required?** 

For the mex compilation, a C compiler must be installed and correctly interfaced to Matlab. Use the command `mex -setup` in Matlab to check your current compiler. The software was successfully tested in Matlab 2016b with the following compilers:

* Linux: GCC compiler
* Mac: clang compiler via XCode 7
* Windows: Microsoft SDK 7.1, Intel Parallel Studio XE 2016, Visual Studio C++ 2015

