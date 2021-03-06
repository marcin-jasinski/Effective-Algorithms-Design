#include "pch.h"
#include "City.h"
#include "BruteForce.h"
#include "BranchAndBound.h"
#include "Stopwatch.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

int citiesNumber;
int** edgesMatrix;

int** readTSPData(std::string);

void runBruteForceBenchmarks();
void runBranchAndBoundBenchmarks();

int main()
{
	edgesMatrix = readTSPData("Test_data/tsp_12.txt");

	// bruteForce(citiesNumber, edgesMatrix);

	branchAndBound(citiesNumber, edgesMatrix);

	// runBruteForceBenchmarks();
	// runBranchAndBoundBenchmarks();
}

// This method reads the matrix data from given data file.
// Returns a two-dimensinal array.
int** readTSPData(std::string fileName)
{
	// Checking given file name format
	if (fileName.find(".txt") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	// Creating the data file handler
	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	citiesNumber = 0;
	int** edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		sourceDataFile >> citiesNumber;
		std::cout << "Number of cities: " << citiesNumber << std::endl;

		// Reading data from the data file based on read number of cities.
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

		std::cout << "\nEdges read from test file: " << std::endl;
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				std::cout << std::setw(4) << edgesMatrix[i][j];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	return edgesMatrix;
}

void runBruteForceBenchmarks()
{
	std::fstream file;
	file.open("Benchmarks/bruteforce_benchmark.txt", std::ios::out);

	if (file.good())
	{
		system("cls");
		file << "\nBrute-force for 6 cities (ver. 1)" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_6_1.txt");
			std::cout << "\nBrute force iteration " << i << " for 6 cities (ver. 1)" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBrute-force for 6 cities (ver. 2)" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_6_2.txt");
			std::cout << "\nBrute force iteration " << i << " for 6 cities (ver. 2)" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBrute-force for 10 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_10.txt");
			std::cout << "\nBrute force iteration " << i << " for 10 cities" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBrute-force for 12 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_12.txt");
			std::cout << "\nBrute force iteration " << i << " for 12 cities" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBrute-force for 13 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_13.txt");
			std::cout << "\nBrute force iteration " << i << " for 13 cities" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBrute-force for 14 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_14.txt");
			std::cout << "\nBrute force iteration " << i << " for 14 cities" << std::endl;

			timer.StartCounter();
			bruteForce(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}
	}
	else
	{
		std::cout << "Error opening benchmark file, abort." << std::endl;
		return;
	}
}

void runBranchAndBoundBenchmarks()
{
	std::fstream file;
	file.open("Benchmarks/branchandbound_benchmark.txt", std::ios::out);

	if (file.good())
	{
		system("cls");
		file << "\nBranch and bound for 6 cities (ver. 1)" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_6_1.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 6 cities (ver. 1)" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBranch and bound for 6 cities (ver. 2)" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_6_2.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 6 cities (ver. 1)" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBranch and bound for 10 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_10.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 10 cities" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBranch and bound for 12 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_12.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 12 cities" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBranch and bound for 13 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_13.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 13 cities" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}

		system("cls");
		file << "\nBranch and bound for 14 cities" << std::endl;
		for (int i = 1; i <= 4; i++)
		{
			double endTime;
			Stopwatch timer = Stopwatch();

			edgesMatrix = readTSPData("Test_data/tsp_14.txt");
			std::cout << "\nBranch and bound iteration " << i << " for 14 cities" << std::endl;

			timer.StartCounter();
			branchAndBound(citiesNumber, edgesMatrix);
			endTime = timer.GetCounter();

			file << endTime << std::endl;
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
		}
	}
	else
	{
		std::cout << "Error opening benchmark file, abort." << std::endl;
		return;
	}
}
