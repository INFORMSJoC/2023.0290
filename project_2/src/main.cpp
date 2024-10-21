
#include <iostream>
#include <array>
#include <ctime>
#include "K_Adapt_Problem_SPP.h"
#include "RobustKSolver.h"

int  main(int argc, char** argv)
{
	int min_delta = 25;
	int max_delta = 75;

	Configuration::ConfigFile* temp_conf = new Configuration::ConfigFile(argv[1]);

	ddid_version = temp_conf->getValueOfKey<bool>("DDID_VERSION");

	delete temp_conf;

	if (!ddid_version)
	{
		min_delta = 100;
		max_delta = 100;
	}

	for (int s = 0; s <= 9; ++s) /*Seed*/
	for (int b = 3; b <= 3; b += 3) /*Budget: 3 and 6*/
	for(int k = 2; k <= 4; ++k) /*Num K*/
			for(int n = 30; n <=50; n +=10) /*Nodes, from 20 to 50 +5*/
					for(int delta = min_delta; delta <=max_delta; delta += 25)
	{

	//int b = 3;// Budget = 3
	// Shortest Path
	SPP data;
	int K = k;
	int B = b;
	int seed = s;
	int N = n;

	//gen_SPP(data, 50, 0,);

	gen_SPP(data, N, seed, K, B, delta);


	K_Adapt_Problem_SPP  info;

	info.setInstanceData(data);
	info.setInstance();

	std::cout<<"	****************	" << std::endl;
	std::cout<<"	Solving instance "<< data.solfilename<<std::endl;
	std::cout<<"	****************	" << std::endl;
	RobustKSolver solver;
	solver.set_configuration_file(argv[1]);
	solver.set_problem_info(&info);
	solver.solve_problem();
	

	info.print_out_file(&solver);

	}

	

	//getchar();
}