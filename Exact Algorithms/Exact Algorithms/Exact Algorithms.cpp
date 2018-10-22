#include "pch.h"
#include "City.h"
#include <iostream>
#include <fstream>

#define EPSILON 0.00000001

City* readSymetricTSPData(std::string);
void readAsymetricTSPData(std::string);

int main()
{
    std::cout << "Hello World!\n"; 
	readSymetricTSPData("TSP_Data/burma14.tsp");
}

City* readSymetricTSPData(std::string fileName)
{
	if (fileName.find(".tsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return NULL;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	// str.erase(0, min(str.find_first_not_of('0'), str.size() - 1));

	int citiesNumber;
	City* citiesArray = new City;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;
		sourceDataFile >> citiesNumber;
		std::cout << "Cities number: " << citiesNumber << std::endl;

		int index;
		double longtitude, latitude;

		while (!sourceDataFile.eof())
		{
			sourceDataFile >> index >> longtitude >> latitude;
			citiesArray[index-1] = City(longtitude + EPSILON, latitude + EPSILON);
		}

		for (int i = 0; i < citiesNumber; i++)
		{
			std::cout << i << " long: " << citiesArray[i].getLongtitude() << " lati: " << citiesArray[i].getLatitude() << " \n";
		}

		return citiesArray;
	}
}

void readAsymetricTSPData(std::string fileName)
{
}
