#ifndef LINK_H
#define LINK_H

#include <Eigen/Dense>
#include <complex>
#include "subgp.h"
using namespace std;
using namespace literals::complex_literals;
using namespace Eigen;

class Link {
	private:
		Matrix3cd A;											// 3x3 matrix
		bool unitary;											// is it a unitary matrix?
	public:
		Link();													// default constructor, identity
		Link(const Matrix3cd initA, const bool isUni=false);	// initializes to a 3x3 matrix
		Link(const SubGpType type, const SubGp sub);			// initialize with 2x2 subgroup
		Link(const Link& U3);									// copy constructor
		Link(const Link&& U3);									// move constructor
		Link inv() const;										// multiplicative inverse
		complex<double> det() const;							// determinant of representation
		complex<double> tr() const;								// trace of representation
		SubGp subGp(const SubGpType type) const;				// get R, S, or T subgp 2x2 matrix
		void gramSchmidt();										// apply Gram Schmidth process
		Link gramSchmidted() const;								// return Gram Schmidth'd link
		void su3Project();										// project to SU3
		Link su3Projected() const;								// return SU3 projected link
		void setA(const Matrix3cd newA, const bool isUni=false);// setter function
		Matrix3cd getA() const;									// getter function
		bool isUnitary() const;
		Link operator=(const Link& U3);
		Link operator=(const Link&& U3);
		Link operator+=(const Link& U3);
		Link operator+(const Link& U3) const;
		Link operator-=(const Link& U3);
		Link operator-(const Link& U3) const;
		Link operator*=(const Link& U3);
		Link operator*(const Link& U3) const;
		Link operator*=(const double& c);
		Link operator*(const double& c) const;
		Link operator*=(const complex<double>& c);
		Link operator*(const complex<double>& c) const;
};

#endif