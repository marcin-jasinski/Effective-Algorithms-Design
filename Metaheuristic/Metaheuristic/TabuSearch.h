#pragma once
#include <vector>

std::vector<int> getRandomPermutationTabu(int);
int pathCost(std::vector<int>, int**, int);

std::vector<int> tabuSearch(int citiesNumber, int** edgesMatrix);
