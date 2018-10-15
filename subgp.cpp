#include "subgp.h"

SubGp::SubGp() {
	A = Matrix2cd::Identity();
	unitary = true;
}

SubGp::SubGp(const Matrix2cd initA, const bool isUni) {
	A = initA;
	unitary = isUni;
}

SubGp::SubGp(const SubGp& U2) {
	A = U2.A;
	unitary = U2.unitary;
}

SubGp::SubGp(const SubGp&& U2) {
	A = U2.A;
	unitary = U2.unitary;
}

SubGp SubGp::inv() const {
	if(unitary) return SubGp(A.adjoint(), true);
	else return SubGp(A.inverse());
}

complex<double> SubGp::det() const {
	return A.determinant();
}

complex<double> SubGp::tr() const {
	return A.trace();
}

void SubGp::gramSchmidt() {
	RowVector2cd u, v;
	u << A(0,0), A(0,1);
	v << A(1,0), A(1,1);

	u.stableNormalize();

	v -= u * u.dot(v);
	v.stableNormalize();
	v -= u * u.dot(v);
	v.stableNormalize();

	A << u,
		 v;

	unitary = true;
}

SubGp SubGp::gramSchmidted() const {
	RowVector2cd u, v;
	u << A(0,0), A(0,1);
	v << A(1,0), A(1,1);

	u.stableNormalize();

	v -= u * u.dot(v);
	v.stableNormalize();
	v -= u * u.dot(v);
	v.stableNormalize();

	Matrix2cd initA;
	initA << u,
			 v;

	return SubGp(initA, true);
}

void SubGp::su2Project() {
	double a0 = (A(0,0) + A(1,1)).real() / 2;
	double a1 = (A(0,1) + A(1,0)).imag() / 2;
	double a2 = (A(0,1) - A(1,0)).real() / 2;
	double a3 = (A(0,0) - A(1,1)).imag() / 2;

	Vector4d proja;
	proja << a0, a1, a2, a3;
	proja.stableNormalize();

	seta(proja);
	unitary = true;
}

SubGp SubGp::su2Projected() const {
	double a0 = (A(0,0) + A(1,1)).real() / 2;
	double a1 = (A(0,1) + A(1,0)).imag() / 2;
	double a2 = (A(0,1) - A(1,0)).real() / 2;
	double a3 = (A(0,0) - A(1,1)).imag() / 2;

	Vector4d proja;
	proja << a0, a1, a2, a3;
	proja.stableNormalize();

	return quaternion(proja);
}

void SubGp::setA(const Matrix2cd newA, const bool isUni) {
	A = newA;
	unitary = isUni;
}

Matrix2cd SubGp::getA() const {
	Matrix2cd copyA = A;
	return copyA;
}

void SubGp::seta(const Vector4d newa) {
	complex<double> a11( newa(0), newa(3));
	complex<double> a12( newa(2), newa(1));
	complex<double> a21(-newa(2), newa(1));
	complex<double> a22( newa(0),-newa(3));

	A << a11, a12,
		 a21, a22;

	unitary = newa.stableNorm() == 1.0;
}

Vector4d SubGp::geta() const {
	double a0 = (A(0,0) + A(1,1)).real() / 2;
	double a1 = (A(0,1) + A(1,0)).imag() / 2;
	double a2 = (A(0,1) - A(1,0)).real() / 2;
	double a3 = (A(0,0) - A(1,1)).imag() / 2;

	Vector4d copya;
	copya << a0, a1, a2, a3;

	return copya;
}

bool SubGp::isUnitary() const {
	return unitary;
}

SubGp SubGp::operator=(const SubGp& U2) {
	if (this != &U2) {
		A = U2.A;
		unitary = U2.unitary;
	}
	return *this;
}

SubGp SubGp::operator=(const SubGp&& U2) {
	if (this != &U2) {
		A = U2.A;
		unitary = U2.unitary;
	}
	return *this;
}

SubGp SubGp::operator+=(const SubGp& U2) {
	A += U2.A;
	unitary = false;
	return *this;
}

SubGp SubGp::operator+(const SubGp& U2) const{
	return SubGp(A + U2.A);
}

SubGp SubGp::operator-=(const SubGp& U2) {
	A -= U2.A;
	unitary = false;
	return *this;
}

SubGp SubGp::operator-(const SubGp& U2) const{
	return SubGp(A - U2.A);
}

SubGp SubGp::operator*=(const SubGp& U2) {
	A *= U2.A;
	unitary &= U2.unitary;
	return *this;
}

SubGp SubGp::operator*(const SubGp& U2) const{
	return SubGp(A * U2.A, unitary && U2.unitary);
}

SubGp SubGp::operator*=(const double& c) {
	A *= c;
	unitary &= (c == 1.0);
	return *this;
}

SubGp SubGp::operator*(const double& c) const{
	return SubGp(A * c, unitary && (c == 1.0));
}

SubGp SubGp::operator*=(const complex<double>& c) {
	A *= c;
	unitary &= (abs(c) == 1.0);
	return *this;
}

SubGp SubGp::operator*(const complex<double>& c) const{
	return SubGp(A * c, unitary && (abs(c) == 1.0));
}

Matrix2cd SubGp::s1() {
	Matrix2cd s1;
	s1 << 0.0, 1.0,
		  1.0, 0.0;
	return s1;
}

Matrix2cd SubGp::s2() {
	Matrix2cd s2;
	s2 << 0.0,-1.i,
		  1.i, 0.0;
	return s2;
}

Matrix2cd SubGp::s3() {
	Matrix2cd s3;
	s3 << 1.0, 0.0,
		  0.0,-1.0;
	return s3;
}

SubGpType SubGp::typeR = R;

SubGpType SubGp::typeS = S;

SubGpType SubGp::typeT = T;

SubGp SubGp::quaternion(Vector4d inita) {
	complex<double> a11( inita(0), inita(3));
	complex<double> a12( inita(2), inita(1));
	complex<double> a21(-inita(2), inita(1));
	complex<double> a22( inita(0),-inita(3));

	Matrix2cd initA;
	initA << a11, a12,
			 a21, a22;

	return SubGp(initA, inita.stableNorm() == 1.0);
}