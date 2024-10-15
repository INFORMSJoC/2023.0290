[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# 2023.0290

# Exact and Approximate Schemes for Robust Optimization Problems with Decision Dependent Information Discovery

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
Exact and Approximate Scheme for Robust Optimization Problems with Decision Dependent Information Discovery (https://doi.org/10.1287/ijoc.2023.0290) by R. Paradiso, A. Georghiou, S. Dabia, D. Tönissen. 
The snapshot is based on 
[this SHA](https://github.com/tkralphs/JoCTemplate/commit/f7f30c63adbcb0811e5a133e1def696b74f3ba15) 
in the development repository. 

**Important: This code is being developed on an on-going basis at 
https://github.com/tkralphs/JoCTemplate. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0290

https://doi.org/10.1287/ijoc.2023.0290.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CacheTest,
  author =        {R. Paradiso, A. Georghiu, S. Dabia, D. Tönissen},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Exact and Approximate Schemes for Robust Optimization Problems with Decision Dependent Information Discovery}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0290.cd},
  url =           {https://github.com/INFORMSJoC/2023.0290},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0290},
}  
```

## Description

This software aims to solve Robust Optimization problems with Decision-Dependent Information Discovery. It implements the exact algorithms and approximate formulations described in [CITE THE PAPER?] to solve instances of the Robust Orientation Problem and the Robust Shortest Path Problem with and without Decision-Dependent Information Discovery.

## Building

The code is written in C++ and requires a C++ compiler that supports the C++17 standard. The code has been tested with the Microsoft C++ (MSVC) compiler (v143 toolset) on Windows. To build the code and generate the executable file, open the project with Visual Studio and build the project. By deafault, the executable file will be named ROPEU.exe.

The code requires IBM ILOG CPLEX Optimization Studio 12.10, in particular the CPLEX C API, to be installed on your system. The code has been tested with CPLEX 12.10. 

## Usage

After building the code, you can run the executable file from the command line. The executable file takes as input a configuration file that specifies several parameters for the algorithm and instance setting. The configuration file (a .cfg file) should be in the following format:

RESULTS_FOLDER				      = Results;
SOLVE_K_ADAPT_OP			      = 1;
SOLVE_EXACT_OP				      = 0;
SOLVE_EXACT					        = 0;
TIME_LIMIT					        = 7200;
DDID_VERSION				        = 1;
MAX_NUM_DISCOVERY_FRACTION  = 0.50;
EXACT_GI_USE_SCENARIOS		  = 0;
EXACT_OP_USE_SCENARIOS      = 0;
EXACT_USE_INFO_CUTS			    = 1;
NUM_K						            = 2;
US_TYPE						          = 2;
US_PARAM					          = 0.10;
ALPHA_SIMMETRY_BREAKING		  = 1;
ALPHA_BOUNDS_REDUCTION		  = 1;
BEST_SCENARIO_CUTS			    = 1;
X_TILDE_CUTS				        = 1;
READ_D_MATRIX				        = 0;
ALRIJNE_COLL_NODE_TIME		  = 360;

where RESULTS_FOLDER can be used to set a folder to store the results. TIME_LIMIT set the CPU Time limit in seconds. DDID_VERSION is used to set the version of the problem to be solved: set it to 1 to solve the Decision-Dependent Information Discovery version of the problem, and 0 otherwise. MAX_NUM_DISCOVERY_FRACTION is used to set the maximum number of discovery that can be used in the solution. In the paper, this is the $\delta$ parameter. The other parameters are problem-specific or algorithm-specific and will be explained in the following sections as needed.
 
### Sensor Placement Orienteering Problem (Section 5.1)

To run the code to solve the Sensor Placement Orienteering Problem and reproduce the results in Section 5.1 of the paper, use the following command:

ROPEU.exe <config_file> -i <instance_file>

where <config_file> is the path to the configuration file, <instance_file> is the path to the instance file. To solve the problem using the exact algorithm described in Section 3, set SOLVE_EXACT_OP = 1 in the configuration file. To solve the K-Adaptability version of the problem, set SOLVE_K_ADAPT_OP = 1 and NUM_K to the desired number of K. US_TYPE is used to set the uncertainty set type. Set it to 2 to use the uncertainty set described in Section 5.1 of the paper. US_PARAM is used to set the parameter of the uncertainty set, i.e., parameter $U$ in the paper. Set this value according to the instance you want to solve as explained in Table 10 to reproduce the results. To use the hybrid decomposition described in Section 3.3. of the paper, set EXACT_OP_USE_SCENARIOS = 1. To use the information cuts described in Section 3.1.1 of the paper, set EXACT_USE_INFO_CUTS = 1, and 0 to use the standard Logic Benders Cuts. If the K-Adaptability version of the problem is solved, the parameters ALPHA_SYMMETRY_BREAKING, ALPHA_BOUNDS_REDUCTION, BEST_SCENARIO_CUTS, and X_TILDE_CUTS are used to activate the symmetry breaking inequality (Section 4.1), the strengthened McCormick (Section 4.2), the Optimistic Inequalities (Section 4.3), and the RLT Inequalities (Section 4.4), respectively.



### Case Study: Alrijne Hospital (Section 6)
To reproduce the results in Section 6 of the paper, run the ...


## Results

Figure 1 in the paper shows the results of the multiplication test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/mult-test.png)

Figure 2 in the paper shows the results of the sum test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/sum-test.png)

## Replicating

To replicate the results in [Cite the paper?] generate the executable file and run the following commands in the terminal.


```
make mult-test
```
or
```
python test.py mult
```
To replicate the results in [Figure 2](results/sum-test), do either

```
make sum-test
```
or
```
python test.py sum
```

## Ongoing Development

This code is being developed on an on-going basis at the author's
[Github site](https://github.com/tkralphs/JoCTemplate).

## Support

For support in using this software, submit an
[issue](https://github.com/tkralphs/JoCTemplate/issues/new).
