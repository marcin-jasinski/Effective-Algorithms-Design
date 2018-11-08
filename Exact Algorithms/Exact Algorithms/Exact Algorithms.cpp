#include "pch.h"
#include "City.h"
#include "BruteForce.h"
#include "BranchAndBound.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <windows.h>

int citiesNumber;
int** edgesMatrix;

int** readTSPData(std::string);

int main()
{
	edgesMatrix = readTSPData("Test_data/tsp_6_1.txt");

	bruteForce(citiesNumber, edgesMatrix);

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (edgesMatrix[i][j] == 0) edgesMatrix[i][j] = INT_MAX;
		}
	}

	branchAndBound(citiesNumber, edgesMatrix);
}

int** readTSPData(std::string fileName)
{
	if (fileName.find(".txt") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	citiesNumber = 0;

	int** edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		sourceDataFile >> citiesNumber;
		std::cout << "Number of cities: " << citiesNumber << std::endl;

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
				std::cout << edgesMatrix[i][j] << "  ";
			}
			std::cout << std::endl;
		}
	}

	return edgesMatrix;
}
