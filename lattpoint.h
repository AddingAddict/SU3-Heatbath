#ifndef LATTPOINT_H
#define LATTPOINT_H

#include <array>
#include <string>
using namespace std;

class LattPoint {
	private:
		array<int,4> x;											// 4-vector representation
		int L;													// lattice size
	public:
		LattPoint(const array<int,4> initX, const int initL);	// constructor
		LattPoint(const LattPoint& x);							// copy constructor
		LattPoint(const LattPoint&& x);							// move constructor
		void setX(const array<int,4> newX);						// setter function
		array<int,4> getX() const;								// getter function
		int getL() const;
		LattPoint scale(const int s);							// scale by integer
		string toString() const;
		LattPoint operator=(const LattPoint& U2);
		LattPoint operator=(const LattPoint&& U2);
		LattPoint operator+=(const LattPoint& U2);
		LattPoint operator+(const LattPoint& U2) const;
		LattPoint operator-=(const LattPoint& U2);
		LattPoint operator-(const LattPoint& U2) const;
		bool operator<(const LattPoint& U2) const;
};

#endif