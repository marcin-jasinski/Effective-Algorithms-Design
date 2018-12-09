#pragma once
#include "pch.h"
#include "Point.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#include <vector>
#include <numeric>
#include <chrono>
#include <ctime>
#include "SimulatedAnnealing.h"
#include "TabuSearch.h"
#include "Stopwatch.h"

static int nint(double);
static double convert_to_geo(double);
int TWOD_geo_distance(const Point a, const Point b);

int** readAsymetricTSPData(std::string);
int** readSymetricTSPData(std::string);

void runSampleSA(int startTemperature, double coolingFactor);
void runSampleTS(int maxIterations);
void printPath(std::vector<int> solution);