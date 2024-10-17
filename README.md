[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact and Approximate Schemes for Robust Optimization Problems with Decision Dependent Information Discovery

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE.txt).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Exact and Approximate Scheme for Robust Optimization Problems with Decision Dependent Information Discovery](https://doi.org/10.1287/ijoc.2023.0290) by R. Paradiso, A. Georghiou, S. Dabia, D. Tönissen. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0290

https://doi.org/10.1287/ijoc.2023.0290.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Paradiso2024,
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

This software aims to solve Robust Optimization problems with Decision-Dependent Information Discovery. It implements the exact algorithms and approximate formulations described in [Exact and Approximate Scheme for Robust Optimization Problems with Decision Dependent Information Discovery](https://doi.org/10.1287/ijoc.2023.0290) to solve instances of the Robust Orientation Problem and the Robust Shortest Path Problem with and without Decision-Dependent Information Discovery.

## Building

The code is written in C++ and requires a C++ compiler that supports the C++17 standard. The code has been tested with the Microsoft C++ (MSVC) compiler (v143 toolset) on Windows. To build the code and generate the executable file, open the project with Visual Studio and build the project. By deafault, the executable file will be named ROPEU.exe.

The code requires IBM ILOG CPLEX Optimization Studio 12.10, in particular the CPLEX C API, to be installed on your system. The code has been tested with CPLEX 12.10. 

## Usage

After building the code, you can run the executable file from the command line. The executable file takes as input a configuration file that specifies several parameters for the algorithm and instance setting. The configuration file (a .cfg file) should be in the following format:

```
RESULTS_FOLDER              = Results;\
SOLVE_K_ADAPT_OP            = 1;\
SOLVE_EXACT_OP              = 0;\
SOLVE_EXACT                 = 0;\
TIME_LIMIT                  = 7200;\
DDID_VERSION                = 1;\
MAX_NUM_DISCOVERY_FRACTION  = 0.50;\
EXACT_GI_USE_SCENARIOS      = 0;\
EXACT_OP_USE_SCENARIOS      = 0;\
EXACT_USE_INFO_CUTS         = 1;\
NUM_K                       = 2;\
US_TYPE                     = 2;\
US_PARAM                    = 0.10;\
ALPHA_SIMMETRY_BREAKING     = 1;\
ALPHA_BOUNDS_REDUCTION      = 1;\
BEST_SCENARIO_CUTS          = 1;\
X_TILDE_CUTS                = 1;\
READ_D_MATRIX               = 0;\
ALRIJNE_COLL_NODE_TIME      = 360;\
```

where
 * `RESULTS_FOLDER` can be used to set a folder to store the results
 * `TIME_LIMIT` set the CPU Time limit in seconds
 * `DDID_VERSION` is used to set the version of the problem to be solved: set it to 1 to solve the Decision-Dependent Information Discovery version of the problem, and 0 otherwise.
 * `MAX_NUM_DISCOVERY_FRACTION` is used to set the maximum number of discovery that can be used in the solution. In the paper, this is the $\delta$ parameter.

The other parameters are problem-specific or algorithm-specific and will be explained in the following sections as needed.
 
### Sensor Placement Orienteering Problem (Section 5.1)

To run the code to solve the Sensor Placement Orienteering Problem and reproduce the results in Section 5.1 of the paper, use the following command:
```
ROPEU.exe <config_file> -i <instance_file>
```
where:
 * `<config_file>` is the path to the configuration file
 * `<instance_file>` is the path to the instance file.

Make sure `SOLVE_EXACT = 0` in the configuration file.

To solve the problem using the exact algorithm described in Section 3, set `SOLVE_EXACT_OP = 1` in the configuration file.

To solve the K-Adaptability version of the problem, set `SOLVE_K_ADAPT_OP = 1` and `NUM_K` to the desired number of `K`.

`US_TYPE` is used to set the uncertainty set type. Set it to 2 to use the uncertainty set described in Section 5.1 of the paper.

`US_PARAM` is used to set the parameter of the uncertainty set, i.e., parameter $U$ in the paper. Set this value according to the instance you want to solve as explained in Table 10 to reproduce the results.

To use the hybrid decomposition described in Section 3.3. of the paper, set `EXACT_OP_USE_SCENARIOS = 1`.

To use the information cuts described in Section 3.1.1 of the paper, set `EXACT_USE_INFO_CUTS = 1`, and 0 to use the standard Logic Benders Cuts.

If the K-Adaptability version of the problem is solved, the parameters `ALPHA_SYMMETRY_BREAKING`, `ALPHA_BOUNDS_REDUCTION`, `BEST_SCENARIO_CUTS`, and `X_TILDE_CUTS` are used to activate the symmetry breaking inequality (Section 4.1), the strengthened McCormick (Section 4.2), the Optimistic Inequalities (Section 4.3), and the RLT Inequalities (Section 4.4), respectively.

### Case Study: Alrijne Hospital (Section 6)

To reproduce the results in Section 6 of the paper, run the executable file as described for the Sensor Placement Orienteering Problem. The instance file should be the files in the "Alrijne_Case_Study_Instances" in the "data" folder, and you should set `READ_D_MATRIX = 1` in the configuration file. The other parameters work as for the Sensor Placement Orienteering Problem. However, to use the uncertainty set $\Xi_2$ in section 6.2.1 set US_TYPE = ?.

### Robust Shortest Path Problem (Appendix E)

To run the code to solve the Robust Shortest Path Problem with the exact algorithm and reproduce the results in Appendix E of the paper, use the following command:

```
ROPEU.exe <config_file> -s <number_of_nodes> <seed_min> <seed_max> <budget_of_uncertainty_set>
```

and make sure `SOLVE_EXACT = 1` in the configuration file and all the other `SOLVE_*` parameters are set to 0. The parameters `<number_of_nodes>`, `<seed_min>`, `<seed_max>`, and `<budget_of_uncertainty_set>` are used to generate the instances: the code generates and solve seed_max - seed_min + 1 networks with number_of_nodes nodes, and the budget of the uncertainty set is set to budget_of_uncertainty_set. To use the information cuts described in Section 3.1.1 of the paper, set `EXACT_USE_INFO_CUTS = 1`, and 0 to use the standard Logic Benders Cuts. To use the hybrid decomposition described in Section 3.3. of the paper, set `EXACT_GI_USE_SCENARIOS = 1`.

### K-Adaptability Robust Shortest Path Problem (Appendix E)

To run the code to solve the K-Adaptability Robust Shortest Path Problem and reproduce the results in Appendix E of the paper, use the following command:

```
K_Adaptability_Solver.exe <config_file>
```

where <config_file> is the path to a configuration file that has the following format: 

```
DDID_VERSION                  = 1;\
USE_ALPHA_SYMMETRY_BREAKING   = 1;\
USE_TIGHTER_MCCORMICK         = 1;\
USE_OPTIMISTIC_CUTS           = 1;\
USE_RLT_CUTS                  = 1;\
```

Set `DDID_VERSION = 1` to solve the Robust Shortest Path Problem with Decision-Dependent Information Discovery, and to 0 otherwise. The other parameters are used to activate the symmetry breaking inequality (Section 4.1), the strengthened McCormick (Section 4.2), the Optimistic Inequalities (Section 4.3), and the RLT Inequalities (Section 4.4), respectively. Set these parameters to 1 to use the corresponding cuts, and 0 otherwise. The code will generate and solve the instances described in Appendix E of the paper. The results will be generated in the folder from which the program is executed.

## Results

Tables 3 - 10 in the paper report the results of the experiments.  

## Replicating

To replicate the results generate the executable file and run it in the terminal using the correct configuration and instance file, as explained above.
