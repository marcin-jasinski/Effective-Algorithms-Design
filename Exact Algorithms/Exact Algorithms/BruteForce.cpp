#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "BruteForce.h"
#include "Stopwatch.h"
#include "Array.h"
#include <cstdlib>
#include <ctime>
#include <windows.h>

int bestTotalCost;
Array bestRoute;

void bruteForce(int citiesNumber, int **edgesMatrix)
{
	Array Numbers_array = Array();
	for (int i = 1; i <= citiesNumber; i++) Numbers_array.pushBack(i);

	double endTime;
	bestTotalCost = INT32_MAX;
	bestRoute = NULL;

	Stopwatch timer = Stopwatch();
	timer.StartCounter();
	permutate(Numbers_array, edgesMatrix, 0);
	endTime = timer.GetCounter();
	std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
	std::cout << " Best route cost: " << bestTotalCost << std::endl;
	return;
}

void permutate(Array numbers, int** edgesMatrix, unsigned int index)
{
	if (index == numbers.getSize())
	{
		calculateTotalRouteCost(numbers, edgesMatrix);
		return;
	}

	for (int i = index; i < numbers.getSize(); i++)
	{
		swapElements(numbers, i, index);
		permutate(numbers, edgesMatrix, index + 1);
		swapElements(numbers, i, index);
	}
}

void swapElements(Array array, int indexA, int indexB)
{
	int x = array.get(indexA);
	array.replaceValueOnIndex(indexA, array.get(indexB));
	array.replaceValueOnIndex(indexB, x);
}

void calculateTotalRouteCost(Array numbers, int** edgesMatrix)
{
	int totalCost = 0;

	int startingCity = numbers.get(0);
	int startingCityIndex = startingCity - 1;

	int lastCity = numbers.get(numbers.getSize() - 1);
	int lastCityIndex = lastCity - 1;

	for (int i = 0; i < numbers.getSize() - 1; i++)
	{
		int currentCity = numbers.get(i) - 1;
		int nextOnRoute = numbers.get(i + 1) - 1;

		totalCost += edgesMatrix[currentCity][nextOnRoute];
	}
	totalCost += edgesMatrix[lastCityIndex][startingCityIndex];

	if (totalCost < bestTotalCost)
	{
		bestTotalCost = totalCost;
		bestRoute = numbers;
	}

	return;
}






