#include "pch.h"
#include "Stopwatch.h"
#include <iostream>

Stopwatch::Stopwatch(){
}

Stopwatch::~Stopwatch(){
}

// This method starts the time counter
void Stopwatch::StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li)) std::cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart);

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

// This method gets and returns current time counter state.
// Returned time will be in seconds.
double Stopwatch::GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}
