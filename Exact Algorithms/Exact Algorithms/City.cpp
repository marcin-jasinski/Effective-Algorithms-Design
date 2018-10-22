#include "pch.h"
#include "City.h"

City::City()
{

}

City::City(double longtitude, double latitude)
{
	this->longtitude = longtitude;
	this->latitude = latitude;
}

City::~City()
{
}

double City::getLongtitude()
{
	return this->longtitude;
}

void City::setLongtitude(double longtitude)
{
	this->longtitude = longtitude;
}

double City::getLatitude()
{
	return this->latitude;
}

void City::setLatitude(double latitude)
{
	this->latitude = latitude;
}
