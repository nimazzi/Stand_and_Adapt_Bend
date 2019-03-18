# investment planning

## Brief Description

This code solves stochastic investment planning problems like

![eq1](https://latex.codecogs.com/gif.latex?%5Cunderset%7B%5Cmathbf%7Bx%7D%20%5Cin%20%5Cmathcal%7BX%7D%7D%7B%5Ctext%7Bmin%7D%7D%20%5C%3B%20f%28%5Cmathbf%7Bx%7D%29%20&plus;%20%5Csum_%7Bi%20%5Cin%20%5Cmathcal%7BI%7D%7D%20%5Cpi_i%20%5Chspace%7B1pt%7D%20g%28x_i%2Cc_i%29%2C)

where ![eq2](https://latex.codecogs.com/gif.latex?f%28%5Cmathbf%7Bx%7D%29) yields the expected total investment and fixed cost, ![eq3](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BI%7D) is the set of decision nodes, ![eq4](https://latex.codecogs.com/gif.latex?%5Cpi_i) is the probability associated to node ![eq5](https://latex.codecogs.com/gif.latex?i). Function ![eq6](https://latex.codecogs.com/gif.latex?g%28x_i%2Cc_i%29) gives the cost gives the cost of operating the system over 5 years, and is formulated as

![eq7](https://latex.codecogs.com/gif.latex?g%28x_i%2Cc_i%29%20%3D%20%5Cunderset%7By_i%20%5Cin%20%5Cmathcal%7BY%7D%7D%7B%5Ctext%7Bmin%7D%7D%5C%7B%20c_i%5E%5Ctop%20%5Chspace%7B-2pt%7D%20C%20y_i%20%5Chspace%7B2pt%7D%20%7C%20%5Chspace%7B2pt%7D%20A%20y_i%20%5Cleq%20B%20x_i%20%5C%7D%2C%20%5Cquad%20%5Cforall%20i%20%5Cin%20%5Cmathcal%7BI%7D)

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
