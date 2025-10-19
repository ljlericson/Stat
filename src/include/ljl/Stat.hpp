#include <iostream>
#include <random>
#include <vector>
#include <tuple>
#include <cmath>

namespace ljl::Stat
{
	class Stat
	{
	public:
		using float64_t = double;
		// this is to make the design more readable
		static inline constexpr uint32_t firstSample = 0; 
	public:

		void startSample()
		{
			m_sampleActive = true;
		}

		void endSample()
		{
			m_sampleData.push_back(std::tuple{
				this->getMean(),
				this->getVar(),
				m_numOfElements
			});
			m_sampleActive = false;
			m_sumX = 0;
			m_sumXsqr = 0;
			m_numOfElements = 0;
		}

		void clear()
		{
			m_sampleData.clear();
		}

		void operator<(int64_t newNumber)
		{
			if(m_sampleActive)
			{
				m_sumX += newNumber;
				m_sumXsqr += newNumber * newNumber;
				m_numOfElements++;
			}
		}

//  |========================================|
//  |============| FUNDEMENTALS |============| 
//  |========================================|

	public:
		float64_t getMean()
		{
			if (!m_sampleActive)
				return -1.0;

			return m_sumX / m_numOfElements;
		}

		float64_t getMean(size_t sample)
		{
			return std::get<0>(m_sampleData.at(sample));
		}

		float64_t getVar()
		{
			if (!m_sampleActive)
				return -1.0;

			float64_t mean = this->getMean();
			return (m_sumXsqr / m_numOfElements) - (mean * mean);
		}

		float64_t getVar(size_t sample)
		{
			return std::get<1>(m_sampleData.at(sample));
		}

		float64_t getUnbiasedEstVar()
		{
			return (1.0 / (m_numOfElements - 1)) * (
				m_sumXsqr - ((m_sumX * m_sumX) / m_numOfElements)
				);
		}

		float64_t getStdDev()
		{
			if (!m_sampleActive)
				return -1.0;

			return sqrt(this->getVar());
		}


//  |========================================|
//	|========| NORMAL APPROXIMATION |========| 
//	|========================================|

	public:
		enum class StdDistTail
		{
			left,
			right
		};

	public:
		float64_t NA_normalApproximationProb(float64_t z, StdDistTail tail, size_t sampleNum)
		{
			float64_t mean = std::get<0>(m_sampleData.at(sampleNum));
			float64_t stdDev = sqrt(std::get<1>(m_sampleData.at(sampleNum)));

			float64_t Z = abs((z - mean)) / stdDev;

			switch (tail)
			{
			case StdDistTail::left:
				if (z > mean)
				{
					return this->p_normalCdf(Z);
				}
				else
				{
					return 1.0f - this->p_normalCdf(Z);
				}
			case StdDistTail::right:
				if (z < mean)
				{
					return this->p_normalCdf(Z);
				}
				else
				{
					return 1.0f - this->p_normalCdf(Z);
				}
				break;
			}

			return -1.0f; // something went wrong if you get this
		}

		float64_t NA_normalApproximationProb(float64_t z, StdDistTail tail)
		{
			float64_t mean = 0;
			float64_t stdDev = 0;
			size_t numElements = 0;
			for (auto [xBar, var, numDataPoints] : m_sampleData)
			{
				mean += xBar;
				stdDev += var;
				numElements++;
			}
			mean /= numElements;
			stdDev = sqrt(stdDev / numElements);


			float64_t Z = abs((z - mean)) / stdDev;

			switch (tail)
			{
			case StdDistTail::left:
				if (z > mean)
				{
					return this->p_normalCdf(Z);
				}
				else
				{
					return 1.0f - this->p_normalCdf(Z);
				}
			case StdDistTail::right:
				if (z < mean)
				{
					return this->p_normalCdf(Z);
				}
				else
				{
					return 1.0f - this->p_normalCdf(Z);
				}
				break;
			}

			return -1.0f; // something went wrong if you get this
		}



//	|========================================|
//	|=========| HYPOTHOSIS TESTING |=========|
//  |========================================|
	
	public:
		enum class HypothTestType
		{
			hasIncreased,
			hasDecreased,
			hasChanged
		};

		// Returns the critical significance level or the probabillity
		// of incorectly rejecting the true null hypothosis or in this
		// case the probabillity of incorectly stating a sample has 
		// increased/decreased/changed from the population
		float64_t HY_getCriticalSignificanLevel(HypothTestType testType, size_t controlSample, size_t testSample)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			float64_t 
				mu    = std::get<0>(m_sampleData.at(controlSample)),
				sigma = sqrt(std::get<1>(m_sampleData.at(controlSample))),

				xBar  = std::get<0>(m_sampleData.at(testSample)),
				s     = sqrt(std::get<1>(m_sampleData.at(testSample)))
			;

			size_t n = std::get<2>(m_sampleData.at(testSample));
			float64_t Z = abs((xBar - mu)) / (s / sqrt(n));

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				if (xBar > mu)
				{
					return 1.0f - this->p_normalCdf(Z);
				}
				else
				{
					return this->p_normalCdf(Z);
				}
					
			case HypothTestType::hasDecreased:
				break;
			case HypothTestType::hasChanged:
				break;
			}

			return -1.0f;
		}





	private:
		float64_t p_normalCdf(float64_t Z)
		{
			return 0.5 * (1.0 + std::erf(Z / std::sqrt(2.0)));
		}

	private:
		bool m_sampleActive = false; 
		float64_t m_sumX = 0;
		float64_t m_sumXsqr = 0;
		size_t m_numOfElements = 0;

		std::vector<std::tuple<float64_t, float64_t, size_t>> m_sampleData;
	};
} 