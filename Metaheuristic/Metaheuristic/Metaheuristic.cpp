#include "pch.h"
#include "commons.h"

const double PI = 3.141592;
const double RRR = 6378.388;

int** edgesMatrix;
int citiesNumber;
int correctSolution;

std::vector<int> solution;

int main()
{
	std::srand(time(NULL));

	edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
	correctSolution = 55000;
	
	std::cout << "Best known solution: " << correctSolution << std::endl;
	std::cout << "Number of cities:    " << citiesNumber << std::endl;

	int startTemperature = 200;
	double coolingFactor = 0.99999;

	int maxIterations = 2000;

	runSampleSA(startTemperature, coolingFactor);

	//runSampleTS(maxIterations);
}

// ==============================================================================

void runSampleSA(int startingTemperature, double clf)
{
	Stopwatch timer = Stopwatch();

	std::cout << "\n=== Simulated Annealing ===" << std::endl;
	std::cout << "\nTemperature: " << startingTemperature << ", cooling factor: " << clf << std::endl;

	timer.StartCounter();

	solution = simulatedAnnealing(startingTemperature, clf, citiesNumber, edgesMatrix);

	double endTime = timer.GetCounter();
	double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
	double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;

	std::cout << "\nElapsed time: " << endTime << " seconds." << std::endl;
	std::cout << "Solution:     " << pathCost << std::endl;
	std::cout << "Error:        " << errorRate << "%" << std::endl;
	std::cout << "Path: " << std::endl;
	printPath(solution);
}

void runSampleTS(int maxIterations)
{
	Stopwatch timer = Stopwatch();

	std::cout << "\n=== Tabu Search ===" << std::endl;
	std::cout << "\nMax iterations: " << maxIterations << std::endl;

	timer.StartCounter();
	solution = tabuSearch(maxIterations, citiesNumber, edgesMatrix);
	double endTime = timer.GetCounter();
	double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
	double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;

	std::cout << "\nElapsed time: " << endTime << " seconds." << std::endl;
	std::cout << "Solution:     " << pathCost << std::endl;
	std::cout << "Error:        " << errorRate << "%" << std::endl;
	std::cout << "Path: " << std::endl;
	printPath(solution);
}

void printPath(std::vector<int> solution)
{
	for (int i = 0; i < solution.size(); i++)
	{
		std::cout << solution.at(i) << " -> ";
	}

	std::cout << solution.at(0) << std::endl;
}

// ==============================================================================
double coolingFactor[5] = { 0.99, 0.999, 0.9999, 0.99999, 0.999999 };

int maxIterationsTab[] = { 1000, 1500, 2000, 2500, 5000, 10000 };
int tabuMultiplierTab[] = { 3, 5, 7, 11 };

void benchmarkTS(std::fstream &file, int** edgesMatrix, int citiesNumber)
{
	Stopwatch timer = Stopwatch();
	file << "\nNumber of cities," << citiesNumber << std::endl;
	for (int i = 0; i < 6; i++)
	{
		int maxIterations = maxIterationsTab[i];

		int bestSolution = INT_MAX;
		double averageTime = 0;
		double averageError = 0;

		system("cls");
		std::cout << "\nNumber of cities: " << citiesNumber << std::endl;
		file << "\nMax iterations," << maxIterations << std::endl;
		std::cout << "Max iterations: " << maxIterations << std::endl;
		std::cout << std::endl;
		file << "Average time,best solution,average error rate" << std::endl;
		for (int i = 1; i <= 5; i++)
		{
			timer.StartCounter();
			solution = tabuSearch(maxIterations, citiesNumber, edgesMatrix);
			double endTime = timer.GetCounter();
			averageTime += endTime;
			double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
			if (pathCost < bestSolution) bestSolution = pathCost;
			double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;
			averageError += averageError;
			std::cout << "Elapsed time: " << endTime << " seconds." << std::endl;
			std::cout << "Solution:     " << pathCost << std::endl;
			std::cout << "Error:        " << errorRate << "%" << std::endl;
		}

		file << averageTime / 3 << "," << bestSolution << "," << averageError / 3 << std::endl;
	}
}

void benchmarkSA(std::fstream &file, int** edgesMatrix, int citiesNumber)
{
	Stopwatch timer = Stopwatch();
	file << "\nNumber of cities," << citiesNumber << std::endl;
	for (int i = 0; i < 5; i++)
	{
		double clf = coolingFactor[i];
		for (int startingTemperature = 200; startingTemperature >= 25; startingTemperature -= 25)
		{
			int bestSolution = INT_MAX;
			double averageTime = 0;
			double averageError = 0;

			system("cls");
			std::cout << "\nNumber of cities: " << citiesNumber << std::endl;
			file << "\nTemperature," << startingTemperature << ",cooling factor," << clf << std::endl;
			std::cout << "Temperature: " << startingTemperature << ", cooling factor: " << clf << std::endl;
			std::cout << std::endl;
			file << "Average time,best solution,average error rate" << std::endl;
			for (int i = 1; i <= 5; i++)
			{
				timer.StartCounter();
				solution = simulatedAnnealing(startingTemperature, clf, citiesNumber, edgesMatrix);
				double endTime = timer.GetCounter();
				averageTime += endTime;
				double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
				if (pathCost < bestSolution) bestSolution = pathCost;
				double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;
				averageError += averageError;
				std::cout << "Elapsed time: " << endTime << " seconds." << std::endl;
				std::cout << "Solution:     " << pathCost << std::endl;
				std::cout << "Error:        " << errorRate << "%" << std::endl;
			}

			file << averageTime / 5 << "," << bestSolution << "," << averageError / 5 << std::endl;
		}
	}
}

void runSimulatedAnnealingBenchmarks(int** edgesMatrix, int citiesNumber, int correctSolution)
{
	std::fstream annealingSA_ATSP;
	annealingSA_ATSP.open("annealingBenchmark_ATSP.txt", std::ios::out | std::ios::trunc);
	if (!annealingSA_ATSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		annealingSA_ATSP << "=== ATSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/br17.atsp");
		correctSolution = 39;
		benchmarkSA(annealingSA_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		benchmarkSA(annealingSA_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");
		correctSolution = 1839;
		benchmarkSA(annealingSA_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/kro124p.atsp");
		correctSolution = 36230;
		benchmarkSA(annealingSA_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		benchmarkSA(annealingSA_ATSP, edgesMatrix, citiesNumber);

		annealingSA_ATSP.close();
	}

	std::fstream annealingSA_TSP;
	annealingSA_TSP.open("annealingBenchmark_TSP.txt", std::ios::out | std::ios::trunc);
	if (!annealingSA_TSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		annealingSA_TSP << "=== TSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
		correctSolution = 3323;
		benchmarkSA(annealingSA_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		benchmarkSA(annealingSA_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		benchmarkSA(annealingSA_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr137.tsp");
		correctSolution = 69853;
		benchmarkSA(annealingSA_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		benchmarkSA(annealingSA_TSP, edgesMatrix, citiesNumber);

		annealingSA_TSP.close();
	}
}

void runTabuSearchBenchmarks(int** edgesMatrix, int citiesNumber, int correctSolution)
{
	std::fstream tabuSearch_ATSP;
	tabuSearch_ATSP.open("tabuSearchBenchmark_ATSP.txt", std::ios::out | std::ios::trunc);
	if (!tabuSearch_ATSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		tabuSearch_ATSP << "=== ATSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/br17.atsp");
		correctSolution = 39;
		benchmarkTS(tabuSearch_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		benchmarkTS(tabuSearch_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");
		correctSolution = 1839;
		benchmarkTS(tabuSearch_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/kro124p.atsp");
		correctSolution = 36230;
		benchmarkTS(tabuSearch_ATSP, edgesMatrix, citiesNumber);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		benchmarkTS(tabuSearch_ATSP, edgesMatrix, citiesNumber);

		tabuSearch_ATSP.close();
	}

	std::fstream tabuSearch_TSP;
	tabuSearch_TSP.open("tabuSearchBenchmark_TSP.txt", std::ios::out | std::ios::trunc);
	if (!tabuSearch_TSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		tabuSearch_TSP << "=== TSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
		correctSolution = 3323;
		benchmarkTS(tabuSearch_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		benchmarkTS(tabuSearch_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		benchmarkTS(tabuSearch_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr137.tsp");
		correctSolution = 69853;
		benchmarkTS(tabuSearch_TSP, edgesMatrix, citiesNumber);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		benchmarkTS(tabuSearch_TSP, edgesMatrix, citiesNumber);

		tabuSearch_TSP.close();
	}
}
// ==============================================================================

int** readAsymetricTSPData(std::string fileName)
{
	if (fileName.find(".atsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		sourceDataFile >> citiesNumber;

		for (int i = 1; i <= 5; i++) sourceDataFile >> readLine;

		edgesMatrix = new int*[citiesNumber];
		for (int i = 0; i < citiesNumber; i++) edgesMatrix[i] = new int[citiesNumber];
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				edgesMatrix[i][j] = -1;
			}
		}

		int edgeCost = 0;
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				sourceDataFile >> edgeCost;
				edgesMatrix[i][j] = edgeCost;
			}
		}
	}

	return edgesMatrix;
}

int** readSymetricTSPData(std::string fileName)
{
	if (fileName.find(".tsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		sourceDataFile >> citiesNumber;

		edgesMatrix = new int*[citiesNumber];
		for (int i = 0; i < citiesNumber; i++) edgesMatrix[i] = new int[citiesNumber];
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++) edgesMatrix[i][j] = -1;
		}

		int currentCityIndex;
		double x, y;

		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("NODE_COORD_SECTION") == std::string::npos);

		std::vector<Point> pointsVector;

		for (int i = 0; i < citiesNumber; i++)
		{
			sourceDataFile >> currentCityIndex >> x >> y;
			pointsVector.push_back(Point(x, y));
		}

		for (int i = 0; i < citiesNumber; i++)
		{
			Point currentCity = pointsVector.at(i);

			for (int j = 0; j < citiesNumber; j++)
			{
				Point nextPoint = pointsVector.at(j);
				if (i == j) edgesMatrix[i][j] = 0;
				else edgesMatrix[i][j] = TWOD_geo_distance(currentCity, nextPoint);
			}
		}
	}

	return edgesMatrix;
}

int nint(double d)
{
	return std::floor(d + 0.5);
}

double convert_to_geo(double x)
{
	int deg = nint(x);
	return PI * (deg + 5.0 * (x - deg) / 3.0) / 180.0;
}

int TWOD_geo_distance(const Point a, const Point b)
{
	Point a_geo(convert_to_geo(a.x), convert_to_geo(a.y));
	Point b_geo(convert_to_geo(b.x), convert_to_geo(b.y));

	double q1 = std::cos(a_geo.y - b_geo.y);
	double q2 = std::cos(a_geo.x - b_geo.x);
	double q3 = std::cos(a_geo.x + b_geo.x);

	return (int)RRR * std::acos(0.5 * ((1.0 + q1)*q2 - (1.0 - q1) * q3)) + 1.0;
}