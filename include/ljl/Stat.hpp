//        __  _ __          _____ __        __ 
//       / / (_) /  _ _    / ___// /_____ _/ /_
//      / / / / /  (_|_)   \__ \/ __/ __ `/ __/
//     / / / / /  _ _     ___/ / /_/ /_/ / /_  
//    /_/_/ /_/  (_|_)   /____/\__/\__,_/\__/  
//     /___/

// (C) LINUS ERICSON 2025
// DOCUMENTATION CAN BE FOUND ON GITHUB PAGE
// https://github.com/ljlericson/Stat

#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <expected>
#include <type_traits>

#if __cplusplus < 201703L
#define STAT_NODISCARD
#else
#define STAT_NODISCARD [[nodiscard]]
#endif

#if __cplusplus < 201103L
#define STAT_NOEXCEPT
#else
#define STAT_NOEXCEPT noexcept
#endif

#if __cplusplus < 202302L
#define STAT_REQUIRE_BASE_SAMPLE
#else
#define STAT_REQUIRE_BASE_SAMPLE requires std::is_base_of_v<__BasicSample, T>
#endif

#ifndef STAT_INLINE
#	if defined(_MSC_VER)
#		define STAT_INLINE __forceinline
#	elif defined(__GNUC__) || defined(__clang__)
#		define STAT_INLINE inline __attribute__((always_inline))
#	else
#		define STAT_INLINE inline
#	endif
#endif

namespace ljl::Stat
{
	namespace // priv stuff
	{
		struct __BasicSample
		{
			virtual double getVar() const STAT_NOEXCEPT = 0;

			virtual double getStdDev() const STAT_NOEXCEPT = 0;

			//virtual std::string getSampleName() const = 0;
		};
	}

	STAT_INLINE double normalCdf(double Z) STAT_NOEXCEPT
	{
		return 0.5 * (1.0 + std::erf(Z / std::sqrt(2.0)));
	}

	STAT_INLINE uint64_t factorial(uint8_t n)
	{
		if (n > 20)
			throw std::runtime_error("N is too large for factorial");

		uint64_t result = 1;
		for (uint8_t i = 1; i <= n; ++i)
		{
			result *= i;
		}
		return result;
	}

	STAT_INLINE uint64_t choose(uint8_t n, uint8_t r)
	{
		return factorial(n) / (factorial(r) * factorial(n - r));
	}

//  |========================================|
//	|=============| EXCEPTIONS |=============| 
//	|========================================|

	class SampleTypeMismatchException : public std::exception
	{
	public:
		~SampleTypeMismatchException()
		{
			if(m_errorCollected && !std::uncaught_exceptions())
				std::cout << "e.what() ---> \nSample types mismatched, i.e discrete sample used in place of continuous\n";
		}

		const char* what() const noexcept override
		{
			m_errorCollected = true;
			return "e.what() ---> \nSample types mismatched, i.e discrete sample used in place of continuous\n";
		}

	private:
		mutable bool m_errorCollected = false;
	};

	class HypothosisTestException : public std::exception
	{
	public:
		HypothosisTestException(const std::string& err)
			: m_err(err)
		{	}

		~HypothosisTestException()
		{
			if (!m_errorCollected && std::uncaught_exceptions())
			{
				if (m_err.empty())
					std::cout << "e.what() ---> \nHypothosis test error\n";
				else
					std::cout << m_err << '\n';
			}
		}

		const char* what() const noexcept override
		{
			if (m_err.empty())
				m_err = "e.what() ---> \nHypothosis test error\n";
			else
				m_err = m_err + '\n';

			m_errorCollected = true;
			return m_err.c_str();
		}

	private:
		mutable bool m_errorCollected = false;
		mutable std::string m_err;
	};

	class InvalidParametreException : public std::exception
	{
	public:
		InvalidParametreException(const std::string& err)
			: m_err(err)
		{	}

		~InvalidParametreException()
		{
			if (!m_errorCollected && std::uncaught_exceptions())
			{
				if (m_err.empty())
					std::cout << "e.what() ---> \nInvalid parametre(s) entered\n";
				else
					std::cout << m_err << '\n';
			}
		}

		const char* what() const noexcept override
		{
			if (m_err.empty())
				m_err = "e.what() ---> \nInvalid parametre(s) entered\n";
			else
				m_err = m_err + '\n';

			m_errorCollected = true;
			return m_err.c_str();
		}

	private:
		mutable bool m_errorCollected = false;
		mutable std::string m_err;
	};

	class InvalidSampleTypeException : public std::exception
	{
	public:
		InvalidSampleTypeException(std::string err)
			: m_err(std::move(err))
		{
		}

		~InvalidSampleTypeException()
		{
			if (!m_errorCollected && std::uncaught_exceptions())
			{
				if (m_err.empty())
					std::cout << "e.what() ---> \nInvalid/unsupported sample type used\n";
				else
					std::cout << m_err << '\n';
			}
		}

		const char* what() const noexcept override
		{
			if (m_err.size() == 0)
			{
				m_err = std::move(std::string{ "e.what() ---> Invalid/unsupported sample type used\n" });
			}

			m_errorCollected = true;
			return m_err.c_str();
		}

	private:
		mutable bool m_errorCollected = false;
		mutable std::string m_err;
	};

//  |========================================|
//	|===========| SAMPLE STRUCTS |===========| 
//	|========================================|

	struct ContinuosSample : public __BasicSample
	{
		double getUnbiasedEstVar() const
		{
			return (1.0 / (numOfElements - 1)) * (
				sumXsqr -
				((sumX * sumX) / numOfElements)
				);
		}

		double getMean() const STAT_NOEXCEPT
		{
			return sumX / numOfElements;
		}

		double getVar() const STAT_NOEXCEPT override
		{
			double mean = this->getMean();
			return (sumXsqr / numOfElements) - (mean * mean);
		}

		double getStdDev() const STAT_NOEXCEPT override
		{
			return std::sqrt(this->getVar());
		}

		void operator<<(double num) STAT_NOEXCEPT
		{
			this->sumX += num;
			this->sumXsqr += num * num;
			this->numOfElements++;
		}

		double sumX = 0;
		double sumXsqr = 0;
		size_t numOfElements = 0;
	};

    STAT_NODISCARD 
    ContinuosSample operator+(const ContinuosSample& lhs, const ContinuosSample& rhs) STAT_NOEXCEPT
    {
        ContinuosSample newSample;
        newSample.numOfElements = lhs.numOfElements + rhs.numOfElements;
        newSample.sumXsqr       = lhs.sumXsqr       + lhs.sumXsqr;
        newSample.sumX          = lhs.sumX          + rhs.sumX;
        
        return newSample;
    }

	struct DiscreteSample : public __BasicSample
	{
		double getUnbiasedEstVar() const
		{
			return (1.0 / (numOfElements - 1)) * (
				sumXsqr -
				((sumX * sumX) / numOfElements)
				);
		}

		double getMean() const STAT_NOEXCEPT
		{
			return sumX / numOfElements;
		}

		double getVar() const STAT_NOEXCEPT override
		{
			double mean = this->getMean();
			return (sumXsqr / numOfElements) - (mean * mean);
		}

		double getStdDev() const STAT_NOEXCEPT override
		{
			return std::sqrt(this->getVar());
		}

		void operator<<(size_t num) STAT_NOEXCEPT
		{
			this->sumX += num;
			this->sumXsqr += num * num;
			this->numOfElements++;
		}

		size_t sumX = 0;
		size_t sumXsqr = 0;
		size_t numOfElements = 0;
	};

	struct BinomialSample : public __BasicSample
	{
		enum struct Case
		{
			success,
			fail
		};

		enum struct Order
		{
			fixed,
			any
		};

		BinomialSample(double p)
			: p(p) { }

		double getProportion() const STAT_NOEXCEPT
		{
			return this->p;
		}

		/*
		IMPORTANT: This function will return an error if the sample size is over 20,
		this is because this function uses factorials and over 20! the number is so
		large that a uint64_t overflows (never thought I would write that)
		@return Returns the chance of the sample randomly occuring given the proportion
		*/
		auto getChanceOfSample(Order order) const -> std::expected<double, std::string>
		{
			// compiler gets annoyed otherwise
			uint8_t n = 0;

			switch (order)
			{
			case Order::any:
				// prevent int overflow due to using uint8
				if (numFails + numSuccsess <= 20)
				{
					n = static_cast<uint8_t>(numFails + numSuccsess);
					return choose(n, static_cast<uint8_t>(numSuccsess)) * (std::pow(p, numSuccsess) * (std::pow((1.0 - p), numFails)));
				}
				else
				{
					return std::unexpected("ERROR: SIZE OF SAMPLE IS TOO LARGE TO ESTIMATE CHANCE OF");
				}

			case Order::fixed:
				return (std::pow(p, numSuccsess) * (std::pow((1.0 - p), numFails)));

			default:
				return -1;
			}
		}

		double getExpected() const STAT_NOEXCEPT
		{
			return (numFails + numSuccsess) * p;
		}

		double getPHat() const STAT_NOEXCEPT
		{
			return static_cast<double>(numSuccsess) / static_cast<double>(numSuccsess + numFails);
		}

		double getVar() const STAT_NOEXCEPT override
		{
			size_t n = numSuccsess + numFails;
			return (p * (1.0 - p)) / n; // variance of proportion
		}

		double getStdDev() const STAT_NOEXCEPT override
		{
			return std::sqrt(this->getVar());
		}

		size_t getNumDataPoints() const STAT_NOEXCEPT
		{
			return numSuccsess + numFails;
		}

		void operator<<(Case _case) STAT_NOEXCEPT
		{
			switch (_case)
			{
			case Case::success:
				numSuccsess++;
				break;
			case Case::fail:
				numFails++;
				break;
			}
		}

		size_t numFails = 0;
		size_t numSuccsess = 0;
		double p = 0.0;
	};

//  |========================================|
//	|========| NORMAL APPROXIMATION |========| 
//	|========================================|

	enum class StdDistTail
	{
		left,
		right
	};

	STAT_INLINE double N_normalApproximationProb(double z, StdDistTail tail, const ContinuosSample& sample) STAT_NOEXCEPT
	{
		double mean = sample.getMean();
		double stdDev = sample.getStdDev();

		double Z = abs((z - mean)) / stdDev;

		switch (tail)
		{
		case StdDistTail::left:
			if (z > mean)
			{
				return normalCdf(Z);
			}
			else
			{
				return 1.0f - normalCdf(Z);
			}
		case StdDistTail::right:
			if (z < mean)
			{
				return normalCdf(Z);
			}
			else
			{
				return 1.0f - normalCdf(Z);
			}
			break;
		}

		return -1.0f; // something went wrong if you get this
	}

	STAT_INLINE double N_normalAproxToBinomial(double z, StdDistTail tail, const BinomialSample& sample) STAT_NOEXCEPT
	{
		double mean = sample.getExpected();
		double stdDev = sample.getStdDev();

		double contCorrection = (tail == StdDistTail::left) ? -0.5 : 0.5;

		double Z = abs(((z + contCorrection) - mean)) / stdDev;

		switch (tail)
		{
		case StdDistTail::left:
			if (z > mean)
			{
				return normalCdf(Z);
			}
			else
			{
				return 1.0f - normalCdf(Z);
			}
		case StdDistTail::right:
			if (z < mean)
			{
				return normalCdf(Z);
			}
			else
			{
				return 1.0f - normalCdf(Z);
			}
			break;
		}

		return -1.0f; // something went wrong if you get this
	}

//	|========================================|
//	|=========| HYPOTHOSIS TESTING |=========|
//  |========================================|
	
	enum class HypothTestType
	{
		hasIncreased,
		hasDecreased,
		hasChanged
	};

	enum class PopVarianceEstimationType
	{
		usePopulation,
		useUnbaisedEstimate
	};


	// Returns the critical significance level or the probabillity
	// of incorectly rejecting the true null hypothosis or in this
	// case the probabillity of incorectly stating a sample has 
	// increased/decreased/changed from the population
	template<typename T> 
		STAT_REQUIRE_BASE_SAMPLE
	STAT_INLINE double HY_getCriticalSignificanLevel(HypothTestType testType, PopVarianceEstimationType estType, const T& controlSample, const T& testSample)
	{
		if constexpr (std::is_same_v<T, ContinuosSample>)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			double
				mu = controlSample.getMean(),
				xBar = testSample.getMean(),

				sigma = 0.0
				;

			switch (estType)
			{
			case PopVarianceEstimationType::usePopulation:
				sigma = controlSample.getStdDev();
				break;
			case PopVarianceEstimationType::useUnbaisedEstimate:
				sigma = std::sqrt(testSample.getUnbiasedEstVar());
				break;
			}

			size_t n = testSample.numOfElements;
			double Z = (xBar - mu) / (sigma / std::sqrt(n));

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				return 1.0 - normalCdf(Z);
			case HypothTestType::hasDecreased:
				return normalCdf(Z);
			case HypothTestType::hasChanged:
				return 2.0 * (1.0 - normalCdf(Z));
			}

			return -1.0f; // prevent compiler warning
		}
		else if constexpr (std::is_same_v<T, DiscreteSample)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			double
				mu = controlSample.getMean(),
				xBar = testSample.getMean(),

				sigma = 0.0
				;

			switch (estType)
			{
			case PopVarianceEstimationType::usePopulation:
				sigma = controlSample.getStdDev();
				break;
			case PopVarianceEstimationType::useUnbaisedEstimate:
				sigma = std::sqrt(testSample.getUnbiasedEstVar());
				break;
			}

			size_t n = testSample.numOfElements;
			double contCorrection = 0.0;
			if(testType != HypothTestType::hasChanged)
				contCorrection = (testType == HypothTestType::hasDecreased) ? -0.5 : 0.5;
			else
			{
				
			}
			double Z = (xBar - mu + contCorrection) / (sigma / std::sqrt(n));

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				return 1.0 - normalCdf(Z);
			case HypothTestType::hasDecreased:
				return normalCdf(Z);
			case HypothTestType::hasChanged:
				return 2.0 * (1.0 - normalCdf(Z));
			}

			return -1.0f; // prevent compiler warning
		}
		else if constexpr (std::is_same_v<T, BinomialSample>)
		{
			// mu   = population mean (assumed)
			// s2   = population sd (assumed)
			// xBar = sample mean
			// s2   = sample sd
			double
				p = controlSample.getProportion(),
				pHat = testSample.getPHat(),

				sigma = 0.0
				;

			switch (estType)
			{
			case PopVarianceEstimationType::usePopulation:
				sigma = controlSample.getStdDev();
				break;
			case PopVarianceEstimationType::useUnbaisedEstimate:
				sigma = std::sqrt((pHat * (1.0 - pHat)) / testSample.getNumDataPoints());
				break;
			}

			if (sigma <= 1e-12)
				throw HypothosisTestException{ "Standard deviation is zero cannot compute Z" };

			size_t n = testSample.getNumDataPoints();
			double Z = (pHat - p) / sigma;

			switch (testType)
			{
			case HypothTestType::hasIncreased:
				return 1.0 - normalCdf(Z);
			case HypothTestType::hasDecreased:
				return normalCdf(Z);
			case HypothTestType::hasChanged:
				return 2.0 * (1.0 - normalCdf(std::abs(Z)));
			}

			return -1.0f; // prevent compiler warning
		}
		else
			throw InvalidSampleTypeException{"Invalid sample type given to hypothosis test, supportted types are Binomial and Continuous"};
	}

	STAT_INLINE bool HY_performHypothTest(HypothTestType testType, 
										  PopVarianceEstimationType estType, 
										  double sigLevel, 
										  const ContinuosSample& controlSample, 
										  const ContinuosSample& testSample)
	{
		// mu   = population mean (assumed)
		// s2   = population sd (assumed)
		// xBar = sample mean
		// s2   = sample sd
		double
			mu = controlSample.getMean(),
			xBar = testSample.getMean(),

			sigma = 0.0
			;

		switch (estType)
		{
		case PopVarianceEstimationType::usePopulation:
			sigma = controlSample.getStdDev();
			break;
		case PopVarianceEstimationType::useUnbaisedEstimate:
			sigma = std::sqrt(testSample.getUnbiasedEstVar());
			break;
		}

		size_t n = testSample.numOfElements;
		double Z = (xBar - mu) / (sigma / std::sqrt(n));

		switch (testType)
		{
		case HypothTestType::hasIncreased:
			return sigLevel > (1.0 - normalCdf(Z));
		case HypothTestType::hasDecreased:
			return sigLevel > (normalCdf(Z));
		case HypothTestType::hasChanged:
			return sigLevel > (2.0 * (1.0 - normalCdf(Z)));
		}

		return false;
	}

	template double ljl::Stat::HY_getCriticalSignificanLevel<ContinuosSample>(HypothTestType, PopVarianceEstimationType, const ContinuosSample&, const ContinuosSample&);
	template double ljl::Stat::HY_getCriticalSignificanLevel<BinomialSample>(HypothTestType, PopVarianceEstimationType, const BinomialSample&, const BinomialSample&);

} 