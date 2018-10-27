#include "pch.h"
#include "City.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#define EPSILON 0.00000001

City* readSymetricTSPData(std::string);
int** readAsymetricTSPData(std::string);

int main()
{
	readAsymetricTSPData("ATSP_Data/br17.atsp");
}

City* readSymetricTSPData(std::string fileName)
{
	if (fileName.find(".tsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return NULL;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	int citiesNumber;
	City* citiesArray = nullptr;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		int citiesNumber;
		int currentCityIndex;
		double latitude;
		double longitude;

		sourceDataFile >> citiesNumber;
		citiesArray = new City[citiesNumber];
		std::cout << "Number of cities: " << citiesNumber << std::endl;
		while (!sourceDataFile.eof())
		{
			sourceDataFile >> currentCityIndex >> latitude >> longitude;
			
			currentCityIndex = currentCityIndex - 1;
			City currentCity = City(latitude, longitude);
			citiesArray[currentCityIndex] = currentCity;
		}
	}

	return citiesArray;
}

int** readAsymetricTSPData(std::string fileName)
{
	if (fileName.find(".atsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	int citiesNumber;
	int** edgesMatrix;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		int citiesNumber;
		int currentCityIndex;
		double latitude;
		double longitude;

		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		std::cout << readLine;
		sourceDataFile >> citiesNumber;
		std::cout << "Number of cities: " << citiesNumber << std::endl;

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
