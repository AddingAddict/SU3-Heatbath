#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>
#include <functional>
#include <map>
#include <Eigen/Dense>
#include "link.h"
#include "subgp.h"
#include "lattpoint.h"
using namespace std;
using namespace Eigen;

struct LattLink {
	array<Link,4> U;
};

const int L = 12;			// grid width
const int L4 = L*L*L*L;		// number of lattice points
const int N = L4*4;			// number of link variables
const int initSteps = 30;	// number of sweeps to thermalize first beta value
const int thermSteps = 10;	// number of sweeps to equilibrate system between betas
const int numMeas = 5;		// number of measurements to meake
const int trashSteps = 2;	// number of sweeps between measurements
const int maxLoop = 2;		// maximum wilson loop size considered

map<LattPoint, LattLink> latt;	// link variables, identified by a point and index in 4D

////////////////////////// random number generators
random_device rd;
mt19937 mtunif(rd());
mt19937 mtgaus(rd());
uniform_real_distribution<double> unifDistr(0,1);
normal_distribution<double> gausDistr(0,1);

auto randUnif = bind(unifDistr, mtunif);
auto randGaus = bind(gausDistr, mtgaus);

////////////////////////// function declarations=
void initialize(const double beta);
void updateLink(const double beta);
Link heatBath(const SubGpType type, const Link W, const double beta);
Vector4d sampleA(const double beta, const double k);
void sweep(const double beta);
double wilsonLoop(const int S);


int main(int argc, char* argv[]) {
	// open data file
	ofstream dataFile("su3.data");
	dataFile << "(beta" << flush;
	for(int i = 1; i <= maxLoop; i++)
		dataFile << ", W(" << i << "x" << i << ")_mean, W(" << i << "x" << i << ")_std" << flush;
	dataFile << ")" << '\n' << '\n' << flush;
	
	// run MC algorithm and calculate correlation length for beta values between
	// 0.1 and 3.0, with a step size of 0.1
	for(int i = 1; i <= 8; i++) {
		double beta = 1 * i;
		cout << " Calculating wilson loops with beta = " << beta << '\n' << flush;

		if(i == 1) {
			// create lattice and decide hot vs cold vs mixed start
			initialize(beta);

			// allow the system to come to thermal equilibrium at the given temperature
			cout << " Performing " << initSteps
				<< " steps to thermalize the system " << flush;
			for(int l = 0; l < initSteps; l++) {
				sweep(beta);
				if((l+1) % 5 == 0)
					cout << "." << flush;
			}
		} else {
			// allow the system to come to thermal equilibrium from the last temperature's equilibrium
			cout << " Performing " << thermSteps
				<< " steps to thermalize the system " << flush;
			for(int l = 0; l < thermSteps; l++) {
				sweep(beta);
				if((l+1) % 5 == 0)
					cout << "." << flush;
			}
		}
	
		// perform production steps
		cout << " Done\n Performing production steps " << flush;
		double loop[maxLoop][numMeas] = {0};
		double loopAv[maxLoop] = {0};
		for(int l = 0; l < numMeas; l++) {
			// sweep 4 times, only measure the last sweep
			for(int n = 0; n < trashSteps; n++)
				sweep(beta);
			
			for(int s = 0; s < maxLoop; s++) {
				loop[s][l] = wilsonLoop(s+1);
				loopAv[s] += loop[s][l];
			}
			cout << "." << flush;
		}
		cout << " Done\n\n" << flush;
		for(int s = 0; s < maxLoop; s++)
			loopAv[s] /= numMeas;
	
		// calculate statistics of correlation length
		double loopStd[maxLoop] = {0};
		for(int s = 0; s < maxLoop; s++) {
			for(int l = 0; l < numMeas; l++) {
				double err = loopAv[s] - loop[s][l];
				loopStd[s] += err*err;
			}

			loopStd[s] /= numMeas - 1;
			loopStd[s] = sqrt(loopStd[s]) / sqrt(numMeas);
		}
		
		// write data to file
		dataFile << "(" << beta << flush;
		for(int s = 0; s < maxLoop; s++)
			dataFile << ", " << loopAv[s] << ", " << loopStd[s] << flush;
		dataFile << ")" << '\n' << flush;
	}
	
	dataFile.close();
}

void initialize(const double beta) {
	// generate a half cold, half hot starting configuration
	for(int i = 0; i < L4; i++) {
		// define Lattpoint
		array<int,4> xi;

		xi[0] =  i        % 10;
		xi[1] = (i/10   ) % 10;
		xi[2] = (i/100  ) % 10;
		xi[3] = (i/1000 ) % 10;

		LattPoint Xi(xi, L);

		// prepare link variable per index
		LattLink Ui;

		for(int mu = 0; mu < 4; mu++) {
			Link UiMu;

			// for a low beta, use a hot start
			if(beta < 4) {
				UiMu.setA(Matrix3cd::Random());
				UiMu.gramSchmidt();
			}
			// for a moderate beta, use a mixed start
			else if(beta < 7) {
				if(randUnif() > 0.5) {
					UiMu.setA(Matrix3cd::Random());
					UiMu.gramSchmidt();
				}
			}
			// otherwise, use a cold start

			Ui.U[mu] = UiMu;
		}

		latt[Xi] = Ui;
	}
}

void updateLink(const double beta, const int i) {
	// pick link variable labeled by integer i
	array<int,4> xi;
	int mu;

	xi[0] =  i        % 10;
	xi[1] = (i/10   ) % 10;
	xi[2] = (i/100  ) % 10;
	xi[3] = (i/1000 ) % 10;
	mu    = (i/10000) % 4;

	// create random lattpoint
	LattPoint Xi(xi, L);

	// create index mu direction
	array<int,4> xmu;

	for(int i = 0; i < 4; i++) {
		if(i == mu) xmu[i] = 1;
		else xmu[i] = 0;
	}

	LattPoint Xmu(xmu, L);

	// Utilde is the sum of 6 products of 3 link variables which interact with the random link
	Link Utilde(Matrix3cd::Zero());

	// calculate Utilde
	for(int nu = 0; nu < 4; nu++) {
		if(mu == nu) continue;

		// create index nu directions
		array<int,4> xnu;

		for(int i = 0; i < 4; i++) {
			if(i == nu) xnu[i] = 1;
			else xnu[i] = 0;
		}

		LattPoint Xnu(xnu, L);

		Utilde += latt[Xi+Xmu].U[nu] * latt[Xi+Xnu].U[mu].inv() * latt[Xi].U[nu].inv() +
			latt[Xi+Xmu-Xnu].U[nu].inv() * latt[Xi-Xnu].U[mu].inv() * latt[Xi-Xnu].U[nu];
	}

	// First perform heat bath for each subgroup
	Link W = latt[Xi].U[mu] * Utilde;
	Link R = heatBath(SubGp::typeR, W, beta);

	W = R * W;
	Link S = heatBath(SubGp::typeS, W, beta);

	W = S * W;
	Link T = heatBath(SubGp::typeT, W, beta);

	// U' = T S R U
	Link UPrime = T * S * R * latt[Xi].U[mu];
	UPrime.gramSchmidt();
	latt[Xi].U[mu] = UPrime;
}

Link heatBath(const SubGpType type, const Link W, const double beta) {
	Vector4d w = W.subGp(type).geta();
	double k = sqrt(abs(W.subGp(type).det()));
	SubGp wbar = SubGp::quaternion(w.stableNormalized());

	SubGp x;

	if(k != 0.0) {
		SubGp xw = SubGp::quaternion(sampleA(beta * 2/3, k));
		x = xw * wbar.inv();
	} else cout << "ah: " << '\n' << flush;

	return Link(type, x);
}

Vector4d sampleA(const double beta, const double k) {
	double a0, a1, a2, a3;

	double w = exp(-2 * beta * k);

	// choose a0 with P(a0) ~ sqrt(1 - a0^2) * exp(beta * k * a0)
	do {
		double xtrial = randUnif()*(1.0 - w) + w;
		a0 = 1 + log(xtrial)/(beta * k);
	} while(sqrt(1-a0*a0) < randUnif());

	// choose \vec{a} randomly
	double r = sqrt(1-a0*a0);

	a1 = randGaus();
	a2 = randGaus();
	a3 = randGaus();

	double norm = sqrt(a1*a1 + a2*a2 + a3*a3);

	a1 *= r / norm;
	a2 *= r / norm;
	a3 *= r / norm;

	Vector4d a;
	a << a0, a1, a2, a3;
	a.stableNormalize();

	return a;
}

void sweep(const double beta){
	// create list of integers labeling each link variable
	vector<int> I(N);
	for(int i = 0; i < N; i++) I[i] = i;

	// shuffle the list
	random_shuffle(I.begin(), I.end());

	// go sequentially through the list, so a sweep hits each 
	for(int i = 0; i < N; i++){
		updateLink(beta, I[i]);
	}
}

double wilsonLoop(const int S) {
	double loopSum = 0;

	// loop over all points in the lattice
	for(int i = 0; i < L4; i++) {
		// define Lattpoint
		array<int,4> xi;

		xi[0] =  i        % 10;
		xi[1] = (i/10   ) % 10;
		xi[2] = (i/100  ) % 10;
		xi[3] = (i/1000 ) % 10;

		LattPoint Xi(xi, L);

		// 6 unique loops in the positive direction at any point
		for(int mu = 0; mu < 4; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				// create index directions
				array<int,4> xmu;
				array<int,4> xnu;

				for(int i = 0; i < 4; i++) {
					if(i == mu) xmu[i] = 1;
					else xmu[i] = 0;

					if(i == nu) xnu[i] = 1;
					else xnu[i] = 0;
				}

				LattPoint Xmu(xmu, L);
				LattPoint Xnu(xnu, L);

				Link loop;

				// build our loop from the 4*S links on the perimeter
				for(int s = 0; s < S; s++)
					loop *= latt[Xi+Xmu.scale(s)].U[mu];
				for(int s = 0; s < S; s++)
					loop *= latt[Xi+Xmu.scale(S)+Xnu.scale(s)].U[nu];
				for(int s = S-1; s >= 0; s--)
					loop *= latt[Xi+Xmu.scale(s)+Xnu.scale(S)].U[mu].inv();
				for(int s = S-1; s >= 0; s--)
					loop *= latt[Xi+Xnu.scale(s)].U[nu].inv();
				loopSum += loop.tr().real() / 3;
			}
		}
	}

	return loopSum / (L4 * 6);
}