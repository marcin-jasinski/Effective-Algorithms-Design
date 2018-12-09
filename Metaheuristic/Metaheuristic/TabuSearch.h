#pragma once
#include <vector>

std::vector<int> getRandomPermutationTabu(int);
int pathCostTabu(std::vector<int>, int**, int);

std::vector<int> tabuSearch(int maxIterations, int citiesNumber, int** edgesMatrix);
