#include "link.h"

Link::Link() {
	A = Matrix3cd::Identity();
	unitary = true;
}

Link::Link(Matrix3cd initA, bool isUni) {
	A = initA;
	unitary = isUni;
}

Link::Link(SubGpType type, SubGp sub) {
	Matrix2cd subA = sub.getA();

	switch(type) {
		case R:
			A << subA(0,0), subA(0,1), 0,
				 subA(1,0), subA(1,1), 0,
				 0,			0,		   1;
			break;

		case S:
			A << subA(0,0), 0, subA(0,1),
				 0,			1, 0,
				 subA(1,0), 0, subA(1,1);
			break;

		case T:
			A << 1, 0,		   0,
				 0, subA(0,0), subA(0,1),
				 0, subA(1,0), subA(1,1);
			break;
	}

	unitary = sub.isUnitary();
}

Link::Link(const Link& U3) {
	A = U3.A;
	unitary = U3.unitary;
}

Link::Link(const Link&& U3) {
	A = U3.A;
	unitary = U3.unitary;
}

Link Link::inv() const {
	if(unitary) return Link(A.adjoint(), true);
	else return Link(A.inverse());
}

complex<double> Link::det() const {
	return A.determinant();
}

complex<double> Link::tr() const {
	return A.trace();
}

SubGp Link::subGp(const SubGpType type) const {
	Matrix2cd subA;

	switch(type) {
		case R:
			subA << A(0,0), A(0,1),
					A(1,0), A(1,1);
			break;
		case S:
			subA << A(0,0), A(0,2),
					A(2,0), A(2,2);
			break;
		case T:
			subA << A(1,1), A(1,2),
					A(2,1), A(2,2);
			break;
	}

	return SubGp(subA, unitary);
}

void Link::gramSchmidt() {
	RowVector3cd u, v, uxv;
	u << A(0,0), A(0,1), A(0,2);
	v << A(1,0), A(1,1), A(1,2);

	u.stableNormalize();

	v -= u * u.dot(v);
	v.stableNormalize();
	v -= u * u.dot(v);
	v.stableNormalize();

	uxv = u.cross(v);
	uxv.stableNormalize();

	A << u,
		 v,
		 uxv;

	unitary = true;
}

Link Link::gramSchmidted() const {
	RowVector3cd u, v, uxv;
	u << A(0,0), A(0,1), A(0,2);
	v << A(1,0), A(1,1), A(1,2);

	u.stableNormalize();

	v -= u * u.dot(v);
	v.stableNormalize();
	v -= u * u.dot(v);
	v.stableNormalize();

	uxv = u.cross(v);
	uxv.stableNormalize();

	Matrix3cd initA;
	initA << u,
			 v,
			 uxv;

	return Link(initA, true);
}

void Link::setA(const Matrix3cd newA, const bool isUni) {
	A = newA;
	unitary = isUni;
}

Matrix3cd Link::getA() const {
	Matrix3cd copyA = A;
	return copyA;
}

bool Link::isUnitary() const {
	return unitary;
}

Link Link::operator=(const Link& U3) {
	if (this != &U3) {
		A = U3.A;
		unitary = U3.unitary;
	}
	return *this;
}

Link Link::operator=(const Link&& U3) {
	if (this != &U3) {
		A = U3.A;
		unitary = U3.unitary;
	}
	return *this;
}

Link Link::operator+=(const Link& U3) {
	A += U3.A;
	unitary = false;
	return *this;
}

Link Link::operator+(const Link& U3) const{
	return Link(A + U3.A);
}

Link Link::operator-=(const Link& U3) {
	A -= U3.A;
	unitary = false;
	return *this;
}

Link Link::operator-(const Link& U3) const{
	return Link(A - U3.A);
}

Link Link::operator*=(const Link& U3) {
	A *= U3.A;
	unitary &= U3.unitary;
	return *this;
}

Link Link::operator*(const Link& U3) const{
	return Link(A * U3.A, unitary && U3.unitary);
}

Link Link::operator*=(const double& c) {
	A *= c;
	unitary &= (c == 1.0);
	return *this;
}

Link Link::operator*(const double& c) const{
	return Link(A * c, unitary && (c == 1.0));
}

Link Link::operator*=(const complex<double>& c) {
	A *= c;
	unitary &= (abs(c) == 1.0);
	return *this;
}

Link Link::operator*(const complex<double>& c) const{
	return Link(A * c, unitary && (abs(c) == 1.0));
}