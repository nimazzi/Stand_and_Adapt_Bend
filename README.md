# Benders with Adaptive Oracles

## Brief Description

This code solves a power system investment planning problem with a time horizon of 15 years. The deterministic version of the problem has 3 decision nodes: one refers to decisions to be taken at present time, one to decisions in 5 years time, and one to decision in 10 years time. The stochastic version is obtained by modeling different possible scenarios for the future of the system in 5 and 10 years. At each node we also compute the cost of operating the system for the following 5 years for given installed capacity. We consider a construction time of 5 years, so new assets installed at present time will only be available in 5 and 10 years, and new capacity installed in 5 years will only be available in 10 years. We model a set ![img1](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BP%7D) of technologies: six thermal units, three storage units, and three renewable generation units.

The stochastic investment planning problem is formulated as

![eq1](https://latex.codecogs.com/gif.latex?%5Cunderset%7B%5Cmathbf%7Bx%7D%20%5Cin%20%5Cmathcal%7BX%7D%7D%7B%5Ctext%7Bmin%7D%7D%20%5C%3B%20f%28%5Cmathbf%7Bx%7D%29%20&plus;%20%5Csum_%7Bi%20%5Cin%20%5Cmathcal%7BI%7D%7D%20%5Cpi_i%20%5Chspace%7B1pt%7D%20g%28x_i%2Cc_i%29%2C)

where ![eq2](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BI%7D) is the set of decision nodes, each associated with a probability ![eq3](https://latex.codecogs.com/gif.latex?%5Cpi_i). The function 

![eq4](https://latex.codecogs.com/gif.latex?g%28x_i%2Cc_i%29%20%3D%20%5Cunderset%7By_i%20%5Cin%20%5Cmathcal%7BY%7D%7D%7B%5Ctext%7Bmin%7D%7D%5C%7B%20c_i%5E%5Ctop%20%5Chspace%7B-2pt%7D%20C%20y_i%20%5Chspace%7B2pt%7D%20%7C%20%5Chspace%7B2pt%7D%20A%20y_i%20%5Cleq%20B%20x_i%20%5C%7D%2C%20%5Cquad%20%5Cforall%20i%20%5Cin%20%5Cmathcal%7BI%7D)

gives the cost of operating the system over 5 years. The vector of right-hand side coefficients ![img2](https://latex.codecogs.com/gif.latex?x_i) is given by 

![eq5](https://latex.codecogs.com/gif.latex?x_i%20%3D%20%5Cleft%28%5Cleft%5C%7Bx%5E%7Bacc%7D_%7Bpi%7D%2C%20%5Cforall%20p%20%5Cin%20%5Cmathcal%7BP%7D%5Cright%5C%7D%2C-%5Cnu%5E%7BD%7D_i%2C%5Cnu%5E%7BE%7D_i%5Cright%29%2C%20%5Cquad%20%5Cforall%20i%20%5Cin%20%5Cmathcal%7BI%7D)

here ![img3](https://latex.codecogs.com/gif.latex?x%5E%7Bacc%7D_%7Bpi%7D) is the accumulated capacity of technology ![img4](https://latex.codecogs.com/gif.latex?p) at node ![img5](https://latex.codecogs.com/gif.latex?i). Parameters ![img6](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BD%7D_i) and ![img7](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i) are the relative level of energy demand and the yearly CO2 emission limit, respectively. The vector of cost coefficients ![img8](https://latex.codecogs.com/gif.latex?c_i) is

![eq6](https://latex.codecogs.com/gif.latex?c_i%20%3D%20%28c%5E%7Bnucl%7D_i%2Cc%5E%7B%5Cmathrm%7Bco%7D_2%7D_i%29%2C%20%5Cquad%20%5Cforall%20i%20%5Cin%20%5Cmathcal%7BI%7D%2C)

where ![img9](https://latex.codecogs.com/gif.latex?c%5E%7Bnucl%7D_i) is the uranium fuel price and ![img10](https://latex.codecogs.com/gif.latex?c%5E%7B%5Cmathrm%7Bco%7D_2%7D_i) the CO2 emission price. Finally, the function ![img11](https://latex.codecogs.com/gif.latex?f%28%5Cmathbf%7Bx%7D%29) yields the expected total investment and fixed cost, and it is computed as

![eq7](https://latex.codecogs.com/gif.latex?f%28%5Cmathbf%7Bx%7D%29%20%3D%20%5Csum_%7Bi%20%5Cin%20%5Cmathcal%7BI%7D%7D%20%5Cpi_%7Bi%7D%20%5Csum_%7Bp%20%5Cin%20%5Cmathcal%7BP%7D%7D%20%5Cleft%28%20c%5E%7Binv%7D_%7Bpi%7D%20x%5E%7Binst%7D_%7Bpi%7D%20&plus;%20c%5E%7Bfix%7D_%7Bpi%7D%20x%5E%7Bacc%7D_%7Bpi%7D%5Cright%29%2C)

where the variable ![img12](https://latex.codecogs.com/gif.latex?x%5E%7Binst%7D_%7Bpi%7D) is the newly installed capacity of technology ![img13](https://latex.codecogs.com/gif.latex?p) at node ![img14](https://latex.codecogs.com/gif.latex?i). Parameters ![img14](https://latex.codecogs.com/gif.latex?c%5E%7Binv%7D_%7Bpi%7D) and ![img15](https://latex.codecogs.com/gif.latex?c%5E%7Bfix%7D_%7Bpi%7D) are the unitary investment and fixed cost of technology ![img16](https://latex.codecogs.com/gif.latex?p) at node ![img14](https://latex.codecogs.com/gif.latex?i). Parameters ![img17](https://latex.codecogs.com/gif.latex?c%5E%7Binv%7D_%7Bpi%7D). The accumulated capacity ![img18](https://latex.codecogs.com/gif.latex?x%5E%7Bacc%7D_%7Bpi%7D) at node ![img19](https://latex.codecogs.com/gif.latex?i) is computed as the sum of the historical capacity and the newly installed capacity ![img21](https://latex.codecogs.com/gif.latex?x%5E%7Binst%7D_%7Bpi%27%7D) in all nodes ![img22](https://latex.codecogs.com/gif.latex?i%27) ancestors to ![img23](https://latex.codecogs.com/gif.latex?i).

We consider three possible sources of uncertainty, i.e., ![img23](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i), ![img24](https://latex.codecogs.com/gif.latex?c%5E%7B%5Cmathrm%7Bco%7D_2%7D_i), and ![img25](https://latex.codecogs.com/gif.latex?c%5E%7Bnucl%7D_i). Each uncertain parameters has 3 possible outcomes in 5 years, each of which is linked to 3 additional possible outcomes in 10 years. The result is 9 possible trajectories for each uncertainty, all with the same probability. We consider 4 different cases of the investment problem. Case 1 is the deterministic version, where ![img26](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i), ![img27](https://latex.codecogs.com/gif.latex?c%5E%7B%5Cmathrm%7Bco%7D_2%7D_i), and ![img28](https://latex.codecogs.com/gif.latex?c%5E%7Bnucl%7D_i) are deterministic parameters (weighted average of the scenarios). Then, case 2 has 1 uncertain parameter, case 3 has 2 uncertain parameters, and case 4 has 3 uncertain parameters. The number of decision nodes for the 4 versions of the problem is

| ![tt1](https://latex.codecogs.com/gif.latex?%5Cmathrm%7Bcase%7D) |  
  ![tt3](https://latex.codecogs.com/gif.latex?%5Cmathrm%7Buncertain%5C%3Bparameters%7D) |
  ![tt2](https://latex.codecogs.com/gif.latex?%7C%5Cmathcal%7BI%7D%7C) | 
|:--:|:--:|:--:|
| 1  | -  | 3  |
| 2  | ![a1](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i) | 13 |
| 3  | ![a2](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i%2Cc%5E%7B%5Cmathrm%7Bco%7D_2%7D_i) | 91 |
| 4  | ![a3](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i%2Cc%5E%7B%5Cmathrm%7Bco%7D_2%7D_i%2Cc%5E%7Bnucl%7D_i)| 757 |

| case | uncertain parameters | decision nodes |
|:----:|:---------------:|:--------------:|
|   0  |        -        |   3            |
|   1  | ![a1](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i)        |   13           |
|   2  | ![a2](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i%2Cc%5E%7B%5Cmathrm%7Bco%7D_2%7D_i)        |   91           |
|   3  | ![a3](https://latex.codecogs.com/gif.latex?%5Cnu%5E%7BE%7D_i%2Cc%5E%7B%5Cmathrm%7Bco%7D_2%7D_i%2Cc%5E%7Bnucl%7D_i)|  757           |

The investment problem can be solved with two algorithms: Algorithm 1 (*Stand_Bend*) and Algorithm 2 (*Adapt_Bend*). Both Algorithm *Stand_Bend* and *Adapt_Bend* are presented in ''Benders Decomposition with Adaptive Oracles for Large Scale Optimization'' (http://www.optimization-online.org/DB_HTML/2019/08/7327.html).

## Prerequisites

Install [Julia 1.0](https://julialang.org/downloads/) (with packages [JuMP.jl 0.18](https://github.com/JuliaOpt/JuMP.jl), [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl), [JLD2.jl](https://github.com/JuliaIO/JLD2.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl)) and [Gurobi 7.5](http://www.gurobi.com/downloads/gurobi-optimizer). 

## Running the code

Open the terminal and change directory to the project folder.
```ShellSession
bash$ cd ~/path_to_folder/Stand_and_Adapt_Bend
```
Start Julia and include "main.jl"
```ShellSession
bash$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.0.2 (2018-11-08)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> include("main.jl")
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
 algorithm 1 -> Stand_Bend (Standard Benders)
 algorithm 2 -> Adapt_Bend (Benders with Adaptive Oracles)
 select Benders-type algorithm:
2
```
The algorithm starts and stops once reached the predifined tolerance, e.g.,
```ShellSession
generating datasets...
precompiling code...
 
decomposition algorithm:
**************************************************
 algorithm          : Adapt_Bend
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
