#pragma once
#include <vector>

std::vector<int> getRandomPermutationTabu(int);
int fitness(std::vector<int>, int**, int);

std::vector<int> tabuSearch(int citiesNumber, int** edgesMatrix);
