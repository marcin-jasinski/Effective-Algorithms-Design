#pragma once
class City
{
private:

	double longtitude;
	double latitude;

public:

	City();
	City(double, double);
	~City();

	double getLongtitude();
	void setLongtitude(double);
	double getLatitude();
	void setLatitude(double);
};

