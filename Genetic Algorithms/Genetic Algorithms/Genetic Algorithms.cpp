#include "pch.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <random>

std::vector<int> getNextSpecimen(int);
std::vector<std::vector<int>> getStartingPopulation(int, int);

int main()
{
	srand(time(NULL));
    std::cout << "Hello World!\n"; 

	std::getchar();
	return 0;
}

std::vector<int> getNextSpecimen(int problemSize)
{
	int genome;
	std::vector<int> specimen;

	for (int i = 1; i <= problemSize; i++)
	{
		do
		{
			genome = std::rand() % problemSize + 1;
		} while (std::find(specimen.begin(), specimen.end(), genome) != specimen.end());

		specimen.push_back(genome);
	}
	
	return specimen;
}

std::vector<std::vector<int>> getStartingPopulation(int populationSize, int problemSize)
{
	std::vector<std::vector<int>> population;

	for (int i = 1; i <= populationSize; i++)
	{
		std::vector<int> newSpecimen = getNextSpecimen(problemSize);
		population.push_back(newSpecimen);
	}

	return population;
}