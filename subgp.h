#ifndef SUBGP_H
#define SUBGP_H

#include <Eigen/Dense>
#include <complex>
using namespace std;
using namespace literals::complex_literals;
using namespace Eigen;

enum SubGpType {R, S, T};

class SubGp {
	private:
		Matrix2cd A;											// 2x2 matrix
		bool unitary;											// is it a unitary matrix?
	public:
		SubGp();												// constructor
		SubGp(const Matrix2cd initA, const bool isUni=false);	// initializes to a 2x2 matrix
		SubGp(const SubGp& U2);									// copy constructor
		SubGp(const SubGp&& U2);								// move constructor
		SubGp inv() const;										// multiplicative inverse
		complex<double> det() const;							// determinant of representation
		complex<double> tr() const;								// trace of representation
		void gramSchmidt();										// apply Gram Schmidth process
		SubGp gramSchmidted() const;							// return Gram Schmidth'd subgp
		void su2Project();										// project to SU2
		SubGp su2Projected() const;								// return SU2 projected subgp
		void setA(const Matrix2cd newA, const bool isUni=false);// setter function (2x2 matrix)
		Matrix2cd getA() const;									// getter function (2x2 matrix)
		void seta(const Vector4d newa);							// setter function (quaternion)
		Vector4d geta() const;									// getter function (projected to quaternion)
		bool isUnitary() const;
		SubGp operator=(const SubGp& U2);
		SubGp operator=(const SubGp&& U2);
		SubGp operator+=(const SubGp& U2);
		SubGp operator+(const SubGp& U2) const;
		SubGp operator-=(const SubGp& U2);
		SubGp operator-(const SubGp& U2) const;
		SubGp operator*=(const SubGp& U2);
		SubGp operator*(const SubGp& U2) const;
		SubGp operator*=(const double& c);
		SubGp operator*(const double& c) const;
		SubGp operator*=(const complex<double>& c);
		SubGp operator*(const complex<double>& c) const;
		
		static Matrix2cd s1();									// static references to Pauli matrices
		static Matrix2cd s2();
		static Matrix2cd s3();

		static SubGpType typeR;									// static references to enum SubGpType
		static SubGpType typeS;
		static SubGpType typeT;

		static SubGp quaternion(Vector4d inita);				// convert from quaternion rep to 2xs matrix
};

#endif