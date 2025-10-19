#include <iostream>
#include "include/ljl/Stat.hpp"


int main()
{
	using namespace ljl::Stat;
	Stat stat{};
	stat.startSample();
	stat < 1;
	stat < 2;
	stat < 3;
	stat < 4;
	stat < 5;
	stat.endSample();

	stat.startSample();
	stat < 2;
	stat < 3;
	stat < 4;
	stat < 5;
	stat < 6;
	stat.endSample();


	std::cout << "Variance: " << stat.getVar(Stat::firstSample + 0)
		      << ", Mean: " << stat.getMean(Stat::firstSample + 0) << '\n';

	std::cout << "Chance is greater than 4: " 
		<< stat.NA_normalApproximationProb(4, Stat::StdDistTail::right) << '\n';

	std::cout << "critical significance level for change:" << '\n'
		<< stat.HY_getCriticalSignificanLevel(Stat::HypothTestType::hasIncreased, 0, 1) << '\n';

	return 0;
}