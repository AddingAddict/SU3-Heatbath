#include "lattpoint.h"
#include <stdexcept>

LattPoint::LattPoint(const array<int,4> initX, const int initL) {
	L = initL;
	for(int i = 0; i < 4; i++)
		x[i] = initX[i] % L;
}

LattPoint::LattPoint(const LattPoint& x2) {
	x = x2.x;
	L = x2.L;
}

LattPoint::LattPoint(const LattPoint&& x2) {
	x = x2.x;
	L = x2.L;
}

void LattPoint::setX(const array<int,4> newX) {
	for(int i = 0; i < 4; i++)
		x[i] = newX[i] % L;
}

array<int,4> LattPoint::getX() const {
	array<int,4> copyX = x;
	return copyX;
}

int LattPoint::getL() const {
	return L;
}

LattPoint LattPoint::scale(const int s) {
	array<int,4> newX;
	for(int i = 0; i < 4; i++)
		newX[i] = (x[i] * s) % L;

	LattPoint scaled(newX, L);
	return scaled;
}

string LattPoint::toString() const {
	return "[" + std::to_string(x[0]) + "," + std::to_string(x[1]) + "," + std::to_string(x[2]) + "," + std::to_string(x[3]) + "]";
}

LattPoint LattPoint::operator=(const LattPoint& x2) {
	if (this != &x2) {
		L = x2.L;
		x = x2.x;
	}
	return *this;
}

LattPoint LattPoint::operator=(const LattPoint&& x2) {
	if (this != &x2)
		L = x2.L;
		x = x2.x;
	return *this;
}

LattPoint LattPoint::operator+=(const LattPoint& x2) {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");

	array<int,4> x2x = x2.getX();

	for(int i = 0; i < 4; i++)
		x[i] = (x[i] + x2x[i]) % L;
	return *this;
}

LattPoint LattPoint::operator+(const LattPoint& x2) const {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");
	
	array<int,4> sumX;

	array<int,4> x2x = x2.getX();

	for(int i = 0; i < 4; i++)
		sumX[i] = (x[i] + x2x[i]) % L;
	LattPoint sum(sumX, L);
	return sum;
}

LattPoint LattPoint::operator-=(const LattPoint& x2) {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");

	array<int,4> x2x = x2.getX();

	for(int i = 0; i < 4; i++) {
		x[i] = x[i] - x2x[i];
		while(x[i] < 0)
			x[i] += L;
		x[i] %= L;
	}
	return *this;
}

LattPoint LattPoint::operator-(const LattPoint& x2) const {
	if(L != x2.getL())
		throw std::invalid_argument("different lattice size");
	
	array<int,4> diffX;

	array<int,4> x2x = x2.getX();

	for(int i = 0; i < 4; i++) {
		diffX[i] = x[i] - x2x[i];
		while(diffX[i] < 0)
			diffX[i] += L;
		diffX[i] %= L;
	}
	LattPoint diff(diffX, L);
	return diff;
}

bool LattPoint::operator<(const LattPoint& x2) const {
	array<int,4> x2x = x2.getX();

	for(int i = 0; i < 4; i++) {
		if(x[i] < x2x[i]) return true;
		else if(x[i] > x2x[i]) return false;
	}
	return false;
}