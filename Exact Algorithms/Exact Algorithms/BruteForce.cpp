#include "pch.h"
#include "BruteForce.h"
#include "Stopwatch.h"
#include "Array.h"

#include <ctime>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <windows.h>

int bestTotalCost;
Array bestRoute;

// This method starts the main brute force algorithm execution.
void bruteForce(int citiesNumber, int **edgesMatrix)
{
	double endTime;
	Stopwatch timer = Stopwatch();
	timer.StartCounter();

	Array Numbers_array = Array();
	for (int i = 1; i <= citiesNumber; i++) Numbers_array.pushBack(i);

	bestTotalCost = INT32_MAX;
	bestRoute = NULL;

	permutate(Numbers_array, edgesMatrix, 0);
	endTime = timer.GetCounter();
	
	std::cout << "=== Brute force ===" << std::endl;
	std::cout << " Best route cost: " << bestTotalCost << std::endl;
	std::cout << " Optimal path: " << bestRoute << std::endl;
	std::cout << " Elapsed time: " << endTime << " seconds." << std::endl;

	return;
}

// This method generates all permutations of a given number set. 
// For each permutation a total cost is calculated.
void permutate(Array numbers, int** edgesMatrix, unsigned int index)
{
	if (index == numbers.getSize())
	{
		// Calculating total cost after generating complete permutation.
		
		//std::cout << numbers << std::endl;
		calculateTotalRouteCost(numbers, edgesMatrix);
		return;
	}

	for (int i = index; i < numbers.getSize(); i++)
	{
		swap(numbers, i, index);
		permutate(numbers, edgesMatrix, index + 1);
		swap(numbers, i, index);
	}
}

// This method is used to swap two elements in numbers array.
void swap(Array array, int indexA, int indexB)
{
	int x = array.get(indexA);
	array.replaceValueOnIndex(indexA, array.get(indexB));
	array.replaceValueOnIndex(indexB, x);
}

// This method calcules and returns total cost of a path created from given vertices.
void calculateTotalRouteCost(Array numbers, int** edgesMatrix)
{
	int totalCost = 0;
	int startingCityIndex = numbers.get(0) - 1;
	int lastCityIndex = numbers.get(numbers.getSize() - 1) - 1;

	for (int i = 0; i < numbers.getSize() - 1; i++)
	{
		int currentCity = numbers.get(i) - 1;
		int nextOnRoute = numbers.get(i + 1) - 1;

		totalCost += edgesMatrix[currentCity][nextOnRoute];
	}

	// Adding cost of returning home from the last vertex
	totalCost += edgesMatrix[lastCityIndex][startingCityIndex];

	// Updating best found cost
	if (totalCost < bestTotalCost)
	{
		bestTotalCost = totalCost;
		bestRoute = NULL;
		for (int i = 0; i < numbers.getSize(); i++) bestRoute.pushBack(numbers.get(i) - 1);
		bestRoute.pushBack(numbers.get(0) - 1);
	}

	return;
}






