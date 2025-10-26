# Stat

Stat is a single header file tool I made for processing and performing calculations on samples of data.

Currently there are binomial and continuous data sample types. These samples can be tested using a Z test to get the critical significance level for whether the test sample has increased, decreased or changed from the control sample.

___

### Example

___

```c++
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
```

___

### Example Application

Takes a control sample of fps and then continuously samples fps to find the significance of how much the mean fps has increased or decreased.

___

```c++
void App::Application::sampleFps()
{
	constexpr float controlDurationMs = 5000.0f;
	constexpr float testDelayMs = 3000.0f; // wait after control ends
	constexpr float testDurationMs = 3000.0f;
	constexpr float sampleIntervalMs = 100.0f;

	static float lastSampleTime = 0.0f;

	switch (m_sampleType)
	{
	case SampleType::control:
		if (!m_sampling) 
		{
			std::cout << "\n--- CAPTURING CONTROL SAMPLE ---\n";
			m_fpsSampleTimer.reset();
			lastSampleTime = 0.0f;
			m_sampling = true;
		}
		else 
		{
			float t = m_fpsSampleTimer.sinceStarted();
			if (t - lastSampleTime >= sampleIntervalMs) 
			{
				m_fpsControlSamle << Util::getFps();
				lastSampleTime = t;
			}

			if (t >= controlDurationMs) 
			{
				m_sampling = false;
				m_sampleType = SampleType::test;
				m_fpsSampleTimer.reset();
				std::cout <<   "=== CONTROL SAMPLE COMPLETE  ===\n"
					<< "MEAN: " << m_fpsControlSamle.getMean() << '\n'
					<< "VAR : " << m_fpsControlSamle.getVar() << '\n';
			}
		}
		break;

	case SampleType::test:
		if (m_reTakeControl)
		{
			m_sampling = false;
			m_reTakeControl = false;
			m_fpsTestSamle = {};
			m_fpsControlSamle = {};
			m_sampleType = SampleType::control;
			break;
		}
		float t = m_fpsSampleTimer.sinceStarted();

		if (!m_sampling && t >= testDelayMs) 
		{
			m_fpsSampleTimer.reset();
			lastSampleTime = 0.0f;
			m_sampling = true;
		}
		else if (m_sampling) 
		{
			float t2 = m_fpsSampleTimer.sinceStarted();
			if (t2 - lastSampleTime >= sampleIntervalMs) 
			{
				m_fpsTestSamle << Util::getFps();
				lastSampleTime = t2;
			}

			if (t2 >= testDurationMs) 
			{
				m_sampling = false;
				m_statResult = ljl::Stat::HY_getCriticalSignificanLevel<ljl::Stat::ContinuosSample>(
					ljl::Stat::HypothTestType::hasChanged,
					ljl::Stat::PopVarianceEstimationType::usePopulation,
					m_fpsControlSamle,
					m_fpsTestSamle);

				m_fpsTestSamle = {};

				// Could also display the stat result for debugging purposes using ImGui
				std::cout << "RESULT OF TEST: " << m_statResult << '\n';
			}
		}
		break;
	}
}
```