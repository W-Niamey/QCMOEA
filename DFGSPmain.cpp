#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "QCMOEA_PM.h"
#include "Individual_PM.h"


int main(int argc, char** argv)
{

	int	instances_number = 405;
	int CPU = 5;
	int Reps = 10;
	vector<vector<Individual_PM>> QCMOEAPMFinalAfterRepParetoSet;
	QCMOEAPMFinalAfterRepParetoSet.clear();
	QCMOEAPMFinalAfterRepParetoSet.resize(instances_number);
	QCMOEA_PM qcmoeapm;
	qcmoeapm.RunEvolution("..\\Benchmark\\", CPU, Reps, QCMOEAPMFinalAfterRepParetoSet, 10, 4, instances_number);
	return 0;


}