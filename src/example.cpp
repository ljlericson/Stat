#include <iostream>
#include "include/ljl/Stat.hpp"


int main()
{
	using namespace ljl::Stat;
	Stat stat{};

	stat.startSample("sample1");
	stat << 1;
	stat << 2;
	stat << 3;
	stat << 4;
	stat << 5;
	stat.endSample();

	stat.startSample("sample2");
	stat << 2;
	stat << 3;
	stat << 4;
	stat << 5;
	stat << 6;
	stat.endSample();

	const Stat::Sample& sample2 = stat.getSample("sample2");
	std::cout << "Number of elements in sample2: " << sample2.numOfElements << '\n';



	std::cout << "Variance1: " << stat.getVar("sample1")
		      << ", Mean1: "   << stat.getMean("sample1") << '\n';

	std::cout << "Variance2: " << stat.getVar("sample2")
			  << ", Mean2: "   << stat.getMean("sample2") << '\n';

	std::cout << "Chance is greater than 4: " << '\n'
		<< stat.N_normalApproximationProb(4, Stat::StdDistTail::right, "sample2") << '\n';

	std::cout << "critical significance level for change:" << '\n'
		<< stat.HY_getCriticalSignificanLevel(Stat::HypothTestType::hasIncreased, "sample1", "sample2") << '\n';

	std::cout << "hypoth test at 10% sig level:" << '\n'
		<< stat.HY_performHypothTest(Stat::HypothTestType::hasIncreased, 0.05, "sample1", "sample2") << '\n';

	return 0;
}