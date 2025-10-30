#include <iostream>
#include <fstream>
#include "include/ljl/Stat.hpp"


void example()
{
	using namespace ljl;
	
	Stat::ContinuosSample sample1{};
	sample1 << 1.0;
	sample1 << 2.0;
	sample1 << 3.0;
	sample1 << 4.0;
	sample1 << 5.0;

	Stat::ContinuosSample sample2{};
	sample1 << 1.0;
	sample1 << 2.0;
	sample1 << 3.0;
	sample1 << 4.0;
	sample1 << 6.0;

	Stat::BinomialSample sampleB1{0.25};
	sampleB1 << Stat::BinomialSample::Case::success;
	sampleB1 << Stat::BinomialSample::Case::fail;
	sampleB1 << Stat::BinomialSample::Case::fail;
	sampleB1 << Stat::BinomialSample::Case::success;

	Stat::BinomialSample sampleB2{ 0.25 };
	sampleB2 << Stat::BinomialSample::Case::success;
	sampleB2 << Stat::BinomialSample::Case::success;
	sampleB2 << Stat::BinomialSample::Case::fail;
	sampleB2 << Stat::BinomialSample::Case::success;

	auto result = sampleB2.getChanceOfSample(Stat::BinomialSample::Order::any);

	if (result)
		std::cout << result.value() << '\n';
	else
		std::cout << result.error() << '\n';

	// other style of err handling, both work
	if(Stat::errIncured())
		std::cout << Stat::getError() << '\n';

	std::cout << Stat::HY_getCriticalSignificanLevel<Stat::BinomialSample>(
		Stat::HypothTestType::hasChanged,
		Stat::PopVarianceEstimationType::usePopulation,
		sampleB1,
		sampleB2
	) << '\n';

	return 0;
}
