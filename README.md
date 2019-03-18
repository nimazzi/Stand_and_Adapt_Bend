# investment planning

## About

We formulate the stochastic investment planning problem as
$$ \alpha = \beta $$

## Prerequisites

Install [Julia 1.0](https://julialang.org/downloads/) (with packages [JuMP 0.18](https://github.com/JuliaOpt/JuMP.jl), [Gurobi 0.5](https://github.com/JuliaOpt/Gurobi.jl), [JLD2 0.1](https://github.com/JuliaIO/JLD2.jl), [CSV 0.4](https://github.com/JuliaData/CSV.jl)) and [Gurobi 7.5](http://www.gurobi.com/downloads/gurobi-optimizer). 

## Running the code

Open the terminal and change directory to the project folder.
```Shell
cd ~/path_to_folder/adaptive-oracles
```
Run main.jl with Julia.
```Shell
julia main.jl
```
You will be asked to select the case study (1, 2, 3 or 4),
```ShellSession
loading packages...
loading functions...

 case 1 -> 0 uncertain parameters
 case 2 -> 1 uncertain parameters
 case 3 -> 2 uncertain parameters
 case 4 -> 3 uncertain parameters
 select case study:
2
```
and the algorithm (1 or 2) used to solve the problem
```ShellSession
 algorithm 1 -> standard
 algorithm 2 -> with adaptive oracles
 select Benders-type algorithm:
2
```
The algorithm starts and stops once reached the predifined tolerance, e.g.,
```ShellSession
generating datasets...
precompiling code...
 
decomposition algorithm:
**************************************************
 algorithm          : Benders with adaptive oracles
 case               : 2
 investment  nodes  : 4
 operational nodes  : 12
--------------------------------------------------
 j =    1, Δ =   99.791 %, t =     2.6 s
 j =    2, Δ =   99.541 %, t =     1.6 s
 ...
 j =   73, Δ =    0.012 %, t =     3.9 s
 j =   74, Δ =    0.009 %, t =     2.4 s

***************************************************************************
***************************************************************************
 
 co2 emission limit : uncertain
 co2 emission cost  : deterministic
 uranium cost       : deterministic
 investment nodes   : 4
 operational nodes  : 12
 optimal objective  : 1.382 x 10^11 £
 
 optimal investments @ 0 years
 ----------------
 tech.         ω1
----------------
 coal           0
 coalccs        0
 OCGT           0
 CCGT          13
 diesel         2
 nuclear        0
 pumpL          0
 pumpH          0
 lithium        0
 onwind        19
 offwind        0
 solar          0
 ----------------
 
 optimal investments @ 5 years
 ------------------------------
 tech.         ω1     ω2     ω3
 ------------------------------
 coal           0      0      0
 coalccs        0      0      2
 OCGT           0      0      0
 CCGT           4      1      0
 diesel         3      2      0
 nuclear        6      9     10
 pumpL          0      0      0
 pumpH          0      0      0
 lithium        0      0      0
 onwind         0      0      0
 offwind        0      0      0
 solar          0      0      0
 ------------------------------
 
 computational results:
 -------------------------------------------------------------------
 ϵ-target (%) | iters     time (s)   RMP (%)   SP (%)   Oracles (%)
 -------------------------------------------------------------------
 1.00         | 42        123         0.05      99.83    0.12
 0.10         | 48        139         0.05      99.82    0.12
 0.01         | 74        258         0.05      99.83    0.12
 -------------------------------------------------------------------

***************************************************************************
***************************************************************************

``` 
