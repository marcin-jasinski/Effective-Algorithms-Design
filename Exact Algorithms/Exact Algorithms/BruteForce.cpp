#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "BruteForce.h"
#include "Array.h"
#include <cstdlib>
#include <ctime>
#include <windows.h>

int bestTotalCost;

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

		// std::cout << "\nCost from " << currentCity + 1 << " to " << nextOnRoute + 1 << ": " << edgesMatrix[currentCity][nextOnRoute];
		totalCost += edgesMatrix[currentCity][nextOnRoute];
	}

	//std::cout << "\nCost from " << lastCity << " to " << startingCity << ": " << edgesMatrix[lastCityIndex][startingCityIndex];
	totalCost += edgesMatrix[lastCityIndex][startingCityIndex];

	/*
	std::cout << "\n Total cost for route ";
	for (int i = 0; i < numbers.getSize(); i++)
	{
		std::cout << numbers.get(i) << " -> ";
	}
	std::cout << numbers[0] << " is " << totalCost << std::endl;*/

	if (totalCost < bestTotalCost)
	{
		bestTotalCost = totalCost;
	}

	return;
}

void swapElements(Array array, int indexA, int indexB)
{
	int x = array.get(indexA);
	array.replaceValueOnIndex(indexA, array.get(indexB));
	array.replaceValueOnIndex(indexB, x);
}

void permutate(Array numbers, int** edgesMatrix, unsigned int index)
{
	if (index == numbers.getSize()) // display vector if you have reached last element 
	{
		calculateTotalRouteCost(numbers, edgesMatrix);
		return;
	}

	for (int i = index; i < numbers.getSize(); i++)
	{
		swapElements(numbers, i, index);	//
		permutate(numbers, edgesMatrix, index + 1);		//	swap elements -> go deeper into the array until you swap last 2 elemets -> rollback 
		swapElements(numbers, i, index);	//
	}
}

LARGE_INTEGER li;
double PCFreq = double(li.QuadPart) / 1000000.0;
__int64 CounterStart = 0;

void StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li)) std::cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart);

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

double GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}

void bruteForce(int citiesNumber, int **edgesMatrix)
{
	// creating and filling number vector 
	Array Numbers_array = Array();
	for (int i = 1; i <= citiesNumber; i++) Numbers_array.pushBack(i);

	double endTime;
	bestTotalCost = INT32_MAX;
	// function call
	std::cout << "\nProcessing, please wait... (if data set is large enough, you have really a lot of time now)" << std::endl;
	StartCounter();
	permutate(Numbers_array, edgesMatrix, 0);
	endTime = GetCounter();
	std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
	std::cout << "\n Best route cost: " << bestTotalCost << std::endl;
	return;
}


