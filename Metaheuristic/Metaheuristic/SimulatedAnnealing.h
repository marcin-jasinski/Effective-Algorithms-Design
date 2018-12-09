#pragma once
#include <vector>

std::vector<int> getRandomPermutation(int);
int getRouteCost(std::vector<int>, int**, int);

std::vector<int> simulatedAnnealing(double maxTemperature, double tempChange, int citiesNumber, int** edgesMatrix);

