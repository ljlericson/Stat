#include <string>
#include <memory>
#include <unordered_map>
#include <tuple>
#include <cmath>

#ifndef STAT_NODISCARD
#define STAT_NODISCARD [[nodiscard]]
#endif

namespace ljl::Stat
{
	class Stat
	{
	public:
		struct Sample
		{
			double sumX = 0;
			double sumXsqr = 0;
			size_t numOfElements = 0;

			std::pair<double, double> meanAndVar;
		};

	public:
		void startSample(const std::string& sampleName)
		{
			// checking if a sample already exists under
			// the name given
			if(!m_sampleActive && !m_samples.contains(sampleName))
			{
				// if it doesn't just make a new sample
				m_activeSample = std::make_unique<Sample>();

				m_activeSampleName = sampleName;
				m_sampleActive = true;
			}
			else if(!m_sampleActive)
			{
				// if the sample already exists we just
				// copy the data and treat it as if it
				// was a new sample
				// Later when .endSample() is called
				// the new sample data will be saved
				m_activeSample = std::make_unique<Sample>(); // make new sample
				*m_activeSample = *m_samples.at(sampleName); // then copy the data

				m_activeSampleName = sampleName;
				m_sampleActive = true;
			}
		}

		void endSample()
		{
			if(m_sampleActive)
			{
				m_activeSample->meanAndVar =
					std::pair{
					this->getMean(),
					this->getVar(),
				};

				if (m_samples.contains(m_activeSampleName))
				{
					m_samples[m_activeSampleName] = std::move(m_activeSample);
				}
				else
				{
					m_samples.insert(std::pair{
						m_activeSampleName,
						std::move(m_activeSample)
					});
				}
				m_activeSampleName.clear();
				m_sampleActive = false;
			}
		}

		void deleteSample(const std::string& sampleName)
		{
			if (m_samples.contains(sampleName))
				m_samples.erase(sampleName);
		}

		STAT_NODISCARD
		const Sample& getSample(const std::string& sampleName)
		{
			if (!m_samples.contains(sampleName))
				std::cout << "ERROR: NO SAMPLE FOUND\n";
			// will crash if no sample exists
			// to prevent underfined behaviour
			// from dangling refs
			return *m_samples.at(sampleName);
		}

		void clear()
		{
			if (m_sampleActive)
				this->endSample();

			m_samples.clear();
			m_activeSampleName.clear();
		}

		void operator<<(double newNumber)
		{
			if(m_sampleActive)
			{
				m_activeSample->sumX += newNumber;
				m_activeSample->sumXsqr += newNumber * newNumber;
				m_activeSample->numOfElements++;
			}
		}

//  |========================================|
//  |============| FUNDEMENTALS |============| 
//  |========================================|

	public:
		double getMean()
		{
			if (!m_sampleActive)
				return -1.0;

			return m_activeSample->sumX / m_activeSample->numOfElements;
		}

		double getMean(const std::string& sample)
		{
			return std::get<0>(m_samples.at(sample)->meanAndVar);
		}

		double getVar()
		{
			if (!m_sampleActive)
				return -1.0;

			double mean = this->getMean();
			return (m_activeSample->sumXsqr / m_activeSample->numOfElements) - (mean * mean);
		}

		double getVar(const std::string& sample)
		{
			return std::get<1>(m_samples.at(sample)->meanAndVar);
		}

		double getUnbiasedEstVar()
		{
			return (1.0 / (m_activeSample->numOfElements - 1)) * (
				m_activeSample->sumXsqr - 
				((m_activeSample->sumX * m_activeSample->sumX) / m_activeSample->numOfElements)
				);
		}

		double getStdDev()
		{
			if (!m_sampleActive)
				return -1.0;

			return std::sqrt(this->getVar());
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
		double N_normalApproximationProb(double z, StdDistTail tail, const std::string& sample)
		{
			double mean = std::get<0>(m_samples.at(sample)->meanAndVar);
			double stdDev = std::sqrt(std::get<1>(m_samples.at(sample)->meanAndVar));

			double Z = abs((z - mean)) / stdDev;

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

		double N_normalApproximationProb(double z, StdDistTail tail)
		{
			double mean = 0;
			double stdDev = 0;
			size_t numElements = 0;
			for (const auto& [key, samplePtr] : m_samples)
			{
				mean += samplePtr->meanAndVar.first;
				stdDev += samplePtr->meanAndVar.second;
				numElements++;
			}
			mean /= numElements;
			stdDev = std::sqrt(stdDev / numElements);


			double Z = std::abs((z - mean)) / stdDev;

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
		double HY_getCriticalSignificanLevel(HypothTestType testType, const std::string& controlSample, const std::string& testSample)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			double 
				mu    = std::get<0>(m_samples.at(controlSample)->meanAndVar),
				sigma = std::sqrt(std::get<1>(m_samples.at(controlSample)->meanAndVar)),

				xBar  = std::get<0>(m_samples.at(testSample)->meanAndVar)
			;

			size_t n = m_samples.at(testSample)->numOfElements;
			double Z = (xBar - mu) / (sigma / std::sqrt(n));

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				return 1.0 - this->p_normalCdf(Z);
			case HypothTestType::hasDecreased:
				return this->p_normalCdf(Z);
			case HypothTestType::hasChanged:
				return 2.0 * (1.0 - this->p_normalCdf(Z));
			}

			return -1.0f;
		}

		bool HY_performHypothTest(HypothTestType testType, double sigLevel, const std::string& controlSample, const std::string& testSample)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			double
				mu = std::get<0>(m_samples.at(controlSample)->meanAndVar),
				sigma = std::sqrt(std::get<1>(m_samples.at(controlSample)->meanAndVar)),

				xBar = std::get<0>(m_samples.at(testSample)->meanAndVar)
				;

			size_t n = m_samples.at(testSample)->numOfElements;
			double Z = (xBar - mu) / (sigma / std::sqrt(n));

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				return sigLevel > (1.0 - this->p_normalCdf(Z));
			case HypothTestType::hasDecreased:
				return sigLevel > (this->p_normalCdf(Z));
			case HypothTestType::hasChanged:
				return sigLevel > (2.0 * (1.0 - this->p_normalCdf(Z)));
			}

			return false;
		}

	private:
		double p_normalCdf(double Z)
		{
			return 0.5 * (1.0 + std::erf(Z / std::sqrt(2.0)));
		}

		double p_calcUnbiasedVar()
		{
			return (1.0 / (m_activeSample->numOfElements - 1)) * (
				m_activeSample->sumXsqr -
				((m_activeSample->sumX * m_activeSample->sumX) / m_activeSample->numOfElements)
				);
		}

	private:
		bool m_sampleActive = false; 
		std::string m_activeSampleName;
		std::unique_ptr<Sample> m_activeSample;

		std::unordered_map<std::string, std::unique_ptr<Sample>> m_samples;
	};
} 