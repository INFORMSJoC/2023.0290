#pragma once
#include "dirent.h"
#include "ROPEUsolver.h"
#include "EXAsolver.h"
#include "ExaDCCG.h"
#include "ProblemDef_SPP.h"


int  main(int argc, char** argv)
{
	DIR* dir = NULL;
	struct dirent* ent = NULL;

	Configuration::ConfigFile* cfg_main = new Configuration::ConfigFile(argv[1]);/*use it as argument*/
	bool solve_k_adapt_op, solve_exact_op, solve_exact;

	solve_k_adapt_op  = cfg_main->getValueOfKey<bool>("SOLVE_K_ADAPT_OP");
	solve_exact_op	  =	cfg_main->getValueOfKey<bool>("SOLVE_EXACT_OP");
	solve_exact		  = cfg_main->getValueOfKey<bool>("SOLVE_EXACT");

	int N, smin, smax, budget = 0;


	int oneInst = 0;
	if ((strcmp(argv[2], "-s") == 0))
	{ 
		N		= atoi(argv[3]);
		smin	= atoi(argv[4]);
		smax	= atoi(argv[5]);
		budget	= atoi(argv[6]);
	}
	else {
		if ((strcmp(argv[2], "-i") == 0))
			oneInst = 1;
		else
			oneInst = 0;

		if (oneInst == 0) {
			dir = opendir(argv[2]);
			if (dir == NULL) {
				printf("Error reading the folder....\n");
				return 0;
				//		return errno;
			}
		}
	}
	

	int go = 1;
	char* file = NULL;
	char buff[1024];
	int exit = 0;
	while (exit != 1)
	{
		if ((strcmp(argv[2], "-s") != 0)){

			if (oneInst == 0 && (ent = readdir(dir)) != NULL) {
				printf("Instance: %s\n", ent->d_name);
				exit = 0;
				snprintf(buff, sizeof(buff), "%s%s%s", argv[2], "\\", ent->d_name);
				printf("%s\n", buff);
				file = buff;
				exit = 0;
			}
			else {
				exit = 1;
			}
		if (oneInst == 0 && exit == 1) {
			break;
		}
		else {
			if (oneInst == 1)
				file = argv[3];
			else {
				if (strcmp(ent->d_name, ".") == 0)
					continue;
				if (strcmp(ent->d_name, "..") == 0)
					continue;
			}
		}
	}


		/*
		DO whatever with file
		*/


		if (solve_exact)
		{


			for (int s = smin; s <= smax; ++s)
			{
			
			// Shortest Path
			SPP data;
			int B = budget;
			int seed = s;

			gen_SPP(data, N, seed, B);


			ProblemDef_SPP  prob;

			prob.setInstanceData(data);
			prob.setConfigFile(argv[1]);
			prob.setInstance();

			//std::string inst_name = data.solfilename.substr(0, data.solfilename.find('.'));
			
			//inst_name.erase(std::remove(inst_name.begin(), inst_name.end(), ".log"), inst_name.end());
			std::string inst_name = data.solfilename;
			std::cout << "\n\nInstance: " << inst_name << "\n" << std::endl;
			
			ExaSolverGI exa_gi_solver;
			exa_gi_solver.set_problem(&prob);

			exa_gi_solver.setConfigFile(argv[1]);

			double opt_val = 0;
			exa_gi_solver.run_exact_method(&opt_val);

			prob.print_out_file(&exa_gi_solver);

			exa_gi_solver.clean_data_structure();

			prob.freeMemory();
		}
			exit = 1;
	}
		// Exa prob
		if (solve_exact_op)
		{
			double opt_val;
			EXAsolver exa_solver;
			exa_solver.read_instance(file, argv[1]);


			exa_solver.run_exact_method(&opt_val);

			exa_solver.print_solution_file();

			exa_solver.print_solution_latex();

			exa_solver.clean_data_structure();
			std::cout << "Exact Value = " << std::to_string(opt_val) << std::endl;
			//getchar();
		}



		// K adapt
		if (solve_k_adapt_op)
		{
			ROPEUsolver solver;
			solver.read_instance(file, argv[1]);

			solver.init_data_structures();


			solver.define_uncertainty_set();

			solver.build_and_solve_model();


			/*Here we want to evaluate for given w*/
			EXAsolver exa_solver;
			double exact_evaluation;
			exa_solver.read_instance(file, argv[1]);
			
			std::vector<bool> w_val;
		

			solver.get_w_values(w_val);

			exa_solver.load_and_init(w_val);

			exa_solver.run_sub_problem_solver(&exact_evaluation);

			solver.exact_w_evaluation = exact_evaluation;

			solver.print_solution_file();

			solver.print_solution_latex();

			solver.print_debug_curr_sol();

			solver.clean_data_structure();
		}

	}
	delete cfg_main;
	return 0;
}