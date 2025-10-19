#include <iostream>
#include "include/ljl/Stat.hpp"


int main()
{
	using namespace ljl::Stat;
	Stat stat{};
	stat.startSample();
	stat < 4.6;
	stat < 6.8;
	stat < 5.2;
	stat < 6.2;
	stat < 5.7;
	stat < 7.1;
	stat < 6.3;
	stat < 5.6;
	stat < 7.0;
	stat < 5.8;
	stat < 6.5;
	stat < 7.2;	
	stat.endSample();

	stat.startSample();
	stat < 4.6;
	stat < 6.8;
	stat < 9.2;
	stat < 6.2;
	stat < 5.7;
	stat < 7.5;
	stat < 6.3;
	stat < 5.6;
	stat < 7.0;
	stat < 5.8;
	stat < 8.0;
	stat < 7.3;
	stat.endSample();


	std::cout << "Variance: " << stat.getVar(Stat::firstSample + 0)
		      << ", Mean: " << stat.getMean(Stat::firstSample + 0) << '\n';

	std::cout << "Chance is greater than 4: " 
		<< stat.NA_normalApproximationProb(4, Stat::StdDistTail::right) << '\n';

	std::cout << "critical significance level for change:" << '\n'
		<< stat.HY_getCriticalSignificanLevel(Stat::HypothTestType::hasIncreased, 0, 1) << '\n';

	return 0;
}