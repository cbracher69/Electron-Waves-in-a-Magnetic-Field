/*****************************************************************************************

  Density and Current Profiles (z-direction) in a Purely Magnetic Field
  (using dimensionless parameters, normalized to free-particle densities and currents)

  Written by:  Christian Bracher
  Version   :  1.03
  Date      :  December 6, 2011

*****************************************************************************************/

/* *** INITIALIZATION ROUTINES *** */

/* Required Libraries */
#include "stdafx.h"
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string>
#include<complex>

using namespace std;

/* Version Info */
	const char *VersionSpec = "1.03";
	const char *VersionDate = "December, 2011";

/* *** GLOBAL VARIABLES *** */

/* Mathematical Constants */
	const double Pi = 4.0 * atan(1.0);
	double PiOver2 = Pi / 2.0;
	complex<double> ii = complex<double>(0.0,1.0);

/* *** Physical Parameters *** */

/* Initial Electron Energy [in units of hbar * omega_L] */
	double Epsilon = 50.00;

/* *** Math Control Parameters *** */

/* Accuracy Goal */
	double Accu = 1.0e-10;

/* *** Output Parameters *** */

/* Observation Window [in m] */
	double xMin = 0.001;
	double xMax = 1.1;

/* Displacement Along Field */
	double zDisp = 3.3;

/* Sampling Number */
  	int xSamples = 1500;

/* The output file name */
  	string OutputFile = "At 50 (profiles).txt";

/* *** CLASSICAL CUTOFF *** */

/* Maximum Number of Cyclotron Orbits */
	long MaxOrbit = 25000;

/* *** QUANTUM CUTOFF *** */

/* Maximum Landau Level */
	long MaxLevel = 10000;

/* *** PROGRAM CONTROL *** */

/* Include Tunneling Trajectories in Semiclassical Approximation? */
	bool IncludeTunneling = true;

/* Integrated Current? */
	bool IntegratedDensity = true;

/* Check for Numerical Inaccuracies? */
	bool CheckNumerics = false;

/* *** MATHEMATICAL ROUTINES *** */

/* The Sign of an Argument */
	int signum(double z)
	{
		int Value = 0;

		if (z > 0) ++Value;
		if (z < 0) --Value;

		return Value;
	}

/* *** Trigonometric Functions */
	inline double cot(double x)
	{
	/* *** Calculate cot(x)
	   *** October 2011, written by CB */

	/* Return Function Value */
		return tan(PiOver2 - x);
	}

/* *** The Airy Functions */
	inline double AiryAi(double x);
	inline double AiryAiPrime(double x);
	inline double AiryBi(double x);
	inline double AiryBiPrime(double x);
	inline complex<double> AiryCi(double x);
	inline complex<double> AiryCiPrime(double x);

/* *** Auxiliary Airy Routines */
	inline double AiryAiSeries(double x);
	inline double AiryAiPrimeSeries(double x);
	inline double AiryBiSeries(double x);
	inline double AiryBiPrimeSeries(double x);
	inline double AiryAiSteed(double x);
	inline double AiryAiPrimeSteed(double x);
	inline complex<double> AiryCiNegAsym(double x);
	inline complex<double> AiryCiPrimeNegAsym(double x);
	inline double AiryAiPosAsym(double x);
	inline double AiryAiPrimePosAsym(double x);
	inline double AiryBiPosAsym(double x);
	inline double AiryBiPrimePosAsym(double x);

inline double AiryAiSeries(double x)
{
/* *** Calculate the series expansion of the Airy function Ai(z)
   *** July 26, 2002, written by CB */

	double c1 = .355028053887817;
	double c2 = .258819403792807;
	double Accu = 1e-20;

	double SumF = 0.0;
	double CF   = 1.0;
	double SumG = 0.0;
	double CG   = x;
	double x3   = x*x*x;
	int    k    = 0;

	do
	{
		SumF += CF;
		SumG += CG;
		k += 3;

		CF *= x3 / (k * (k - 1));
		CG *= x3 / (k * (k + 1));
	}
	while ((fabs(CF)+fabs(CG)) > Accu);

	return c1 * SumF - c2 * SumG;
}

inline double AiryBiSeries(double x)
{
/* *** Calculate the series expansion of the Airy function Bi(z)
   *** July 26, 2002, written by CB */

	double c1 = .355028053887817;
	double c2 = .258819403792807;
	double Accu = 1e-20;

	double SumF = 0.0;
	double CF   = 1.0;
	double SumG = 0.0;
	double CG   = x;
	double x3   = x*x*x;
	int    k    = 0;

	do
	{
		SumF += CF;
		SumG += CG;
		k += 3;

		CF *= x3 / (k * (k - 1));
		CG *= x3 / (k * (k + 1));
	}
	while ((fabs(CF)+fabs(CG)) > Accu);

	return sqrt(3.0) * (c1 * SumF + c2 * SumG);
}

inline double AiryAiPrimeSeries(double x)
{
/* *** Calculate the series expansion of the Airy function Ai'(z)
   *** July 26, 2002, written by CB */

	double c1 = .355028053887817;
	double c2 = .258819403792807;
	double Accu = 1e-20;

	double SumFp = 0.0;
	double CFp   = .5*x*x;
	double SumGp = 0.0;
	double CGp   = 1.0;
	double x3    = x*x*x;
	int    k     = 0;

	do
	{
		SumFp += CFp;
		SumGp += CGp;
		k += 3;

		CFp *= x3 / (k * (k + 2));
		CGp *= x3 / (k * (k - 2));
	}
	while ((fabs(CFp) + fabs(CGp)) > Accu);

	return c1 * SumFp - c2 * SumGp;
}

inline double AiryBiPrimeSeries(double x)
{
/* *** Calculate the series expansion of the Airy function Bi'(z)
   *** July 26, 2002, written by CB */

	double c1 = .355028053887817;
	double c2 = .258819403792807;
	double Accu = 1e-20;

	double SumFp = 0.0;
	double CFp   = .5*x*x;
	double SumGp = 0.0;
	double CGp   = 1.0;
	double x3    = x*x*x;
	int    k     = 0;

	do
	{
		SumFp += CFp;
		SumGp += CGp;
		k += 3;

		CFp *= x3 / (k * (k + 2));
		CGp *= x3 / (k * (k - 2));
	}
	while ((fabs(CFp)+fabs(CGp)) > Accu);

	return sqrt(3.0) * (c1 * SumFp + c2 * SumGp);
}

inline double AiryAiSteed(double x)
{
/* *** Calculate the Airy function Ai(z)
   *** for positive z via Steed's algorithm (logarithmic derivative, Wronskian)
   *** July 26, 2002, written by CB */

   double z =.75 / (x * sqrt(x));
   double Accu = 1e-18;

   /* *** Steed's algorithm to calculate Ai'(z)/Ai(z) */
   double Sum = 1.0 + z / 3.0;
   double ak = -5.0 / 144.0;
   double bk = 1.0 + .5 / z;

   double D = 1.0 / bk;
   double Delta = 4.0 * z * ak * D;
   int    k = 1;

   while (fabs(Delta / Sum) > Accu)
   {
	   Sum += Delta;
	   ak -= k / 2.0;
		bk += 1.0;
	   ++k;
		D = 1.0 / (D * ak + bk);
		Delta *= (bk * D - 1.0);
   }

/* Calculate Ai(z) by Wronskian */
	return 1.0 / (Pi * (AiryBiPrimeSeries(x) + sqrt(x) * Sum * AiryBiSeries(x)));
}

inline double AiryAiPrimeSteed(double x)
{
/* *** Calculate the Airy function Ai'(z)
   *** for positive z via Steed's algorithm (logarithmic derivative, Wronskian)
   *** July 26, 2002, written by CB */

   double z =.75 / (x * sqrt(x));
   double Accu = 1e-18;

   /* *** Steed's algorithm to calculate Ai'(z)/Ai(z) */
   double Sum = 1.0 + z / 3.0;
   double ak = -5.0 / 144.0;
   double bk = 1.0 + .5 / z;

   double D = 1.0 / bk;
   double Delta = 4.0 * z * ak * D;
   int    k = 1;

   while (fabs(Delta / Sum) > Accu)
   {
	   Sum += Delta;
	   ak -= k / 2.0;
		++bk;
	   ++k;
		D = 1.0 / (D * ak + bk);
		Delta *= (bk * D - 1.0);
   }

/* Calculate Ai'(z) by Wronskian */
	return -sqrt(x) * Sum / (Pi * (AiryBiPrimeSeries(x) + sqrt(x) * Sum * AiryBiSeries(x)));
}

inline complex<double> AiryCiNegAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Ci(z)
   *** for large negative arguments z
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

   complex<double> iz = ii * 2.0 * (-x) * sqrt(-x) / 3.0;
   complex<double> izInverse = 1.0 / iz;
   complex<double> Sum = 0.0;
   complex<double> Term = 1.0;
   complex<double> Old = 0.0;
   double Aux = 5.0 / 72.0;
   double Accu = 1e-36;
   int k = 0;

/* Calculate series */
   do
   {
		Sum += Term;
		Old = Term;
		Term *= (0.5 * k + Aux / (k + 1)) * izInverse;
		++k;
   }
   while ((norm(Term) < norm(Old)) && (norm(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor */
   Sum *= exp(iz + .25 * ii * Pi) / sqrt(Pi * sqrt(-x));

   return Sum;
}

inline complex<double> AiryCiPrimeNegAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Ci'(z)
   *** for large negative arguments z
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

   complex<double> iz = ii * 2.0 * (-x) * sqrt(-x) / 3.0;
   complex<double> izInverse = 1.0 / iz;
   complex<double> Sum = 0.0;
   complex<double> Term = 1.0;
   complex<double> Old = 0.0;
   double Aux = - 7.0 / 72.0;
   double Accu = 1e-36;
   int k = 0;

/* Calculate series */
   do
   {
		Sum += Term;
		Old = Term;
		Term *= (0.5 * k + Aux / (k + 1)) * izInverse;
		++k;
   }
   while ((norm(Term) < norm(Old)) && (norm(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor and exponential */
	Sum *= (-ii) * sqrt(sqrt(-x) / Pi) * exp(iz + .25 * ii * Pi);

	return Sum;
}

inline double AiryAiPosAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Ai(z)
   *** for large positive arguments z, exponential prefactor included
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

	double zeta = 2.0 * x * sqrt(x) / 3.0;
	double zetaInverse = 1.0 / zeta;
	double Sum = 0.0;
	double Term = 1.0;
	double Old = 0.0;
	double Accu = 1e-18;
	double Aux = 5.0 / 72.0;
	int k = 0;

/* Calculate series */
   do
   {
	   Sum += Term;
	   Old = Term;
	   Term *= -(0.5 * k + Aux / (k + 1)) * zetaInverse;
	   ++k;
   }
   while ((fabs(Term) < fabs(Old)) && (fabs(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor and exponential */
   return .5 * exp(-zeta) * Sum / sqrt(Pi * sqrt(x));
}

inline double AiryAiPrimePosAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Ai'(z)
   *** for large positive arguments z, exponential prefactor included
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

   double zeta = 2.0 * x * sqrt(x) / 3.0;
   double zetaInverse = 1.0 / zeta;
   double Sum = 0.0;
   double Term = 1.0;
   double Old = 0.0;
   double Aux = -7.0 / 72.0;
   double Accu = 1e-18;
   int k = 0;

/* Calculate series */
   do
   {
	   Sum += Term;
	   Old = Term;
	   Term *= -(0.5 * k + Aux / (k + 1)) * zetaInverse;
	   ++k;
   }
   while ((fabs(Term) < fabs(Old)) && (fabs(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor and exponential */
   return -.5 * sqrt(sqrt(x) / Pi) * exp(-zeta) * Sum;
}

inline double AiryBiPosAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Bi(z)
   *** for large positive arguments z, exponential prefactor included
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

   double zeta = 2.0 * x * sqrt(x) / 3.0;
   double zetaInverse = 1.0 / zeta;
   double Sum = 0.0;
   double Term = 1.0;
   double Old = 0.0;
   double Accu = 1e-18;
   double Aux = 5.0 / 72.0;
   int k = 0;

/* Calculate series */
   do
   {
	   Sum += Term;
	   Old = Term;
	   Term *= (0.5 * k + Aux / (k + 1)) * zetaInverse;
	   ++k;
   }
   while ((fabs(Term) < fabs(Old)) && (fabs(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor (not the exponential) */
   return exp(zeta) * Sum / sqrt(Pi * sqrt(x));
}

inline double AiryBiPrimePosAsym(double x)
{
/* *** Calculate the asymptotic series expansion of the Airy function Bi'(z)
   *** for large positive arguments z, exponential prefactor included
   *** July 26, 2002, written by CB; simplified Oct 5, 2011 */

   double zeta = 2.0 * x * sqrt(x) / 3.0;
   double zetaInverse = 1.0 / zeta;
   double Sum = 0.0;
   double Term = 1.0;
   double Old = 0.0;
   double Aux = -7.0 / 72.0;
   double Accu = 1e-18;
   int k = 0;

/* Calculate series */
   do
   {
	   Sum += Term;
	   Old = Term;
	   Term *= (0.5 * k + Aux / (k + 1)) * zetaInverse;
	   ++k;
   }
   while ((fabs(Term) < fabs(Old)) && (fabs(Term) > Accu));

   Sum += .5 * Term;

/* Multiply with prefactor (not the exponential) */
   return exp(zeta) * sqrt(sqrt(x) / Pi) * Sum;
}

inline double AiryAi(double x)
{
/* *** Calculate the Airy function Ai(z)
   *** July 26, 2002, written by CB */

   double Ai;
   if      (x < - 6.2)  Ai = imag(AiryCiNegAsym(x));
   else if (x <   2.0)  Ai = AiryAiSeries(x);
   else if (x <  10.0)  Ai = AiryAiSteed(x);
   else                 Ai = AiryAiPosAsym(x);

   return Ai;
}

inline double AiryAiPrime(double x)
{
/* *** Calculate the Airy function Ai'(z)
   *** July 26, 2002, written by CB */

   double AiP;
   if      (x < - 6.2)  AiP = imag(AiryCiPrimeNegAsym(x));
   else if (x <   2.0)  AiP = AiryAiPrimeSeries(x);
   else if (x <  10.0)  AiP = AiryAiPrimeSteed(x);
   else                 AiP = AiryAiPrimePosAsym(x);

   return AiP;
}

inline double AiryBi(double x)
{
/* *** Calculate the Airy function Bi(z)
   *** July 26, 2002, written by CB */

   double Bi;
   if      (x < - 6.2)  Bi = real(AiryCiNegAsym(x));
   else if (x <  10.0)  Bi = AiryBiSeries(x);
   else                 Bi = AiryBiPosAsym(x);

   return Bi;
}

inline double AiryBiPrime(double x)
{
/* *** Calculate the Airy function Bi'(z)
   *** July 26, 2002, written by CB */

   double BiP;
   if      (x < - 6.2)  BiP = real(AiryCiPrimeNegAsym(x));
   else if (x <  10.0)  BiP = AiryBiPrimeSeries(x);
   else                 BiP = AiryBiPrimePosAsym(x);

   return BiP;
}

inline complex<double> AiryCi(double x)
{
/* *** Calculate the Airy function Ci(z)
   *** July 26, 2002, written by CB */

   complex<double> Ci;
   if      (x < - 6.2)  Ci = AiryCiNegAsym(x);
   else if (x <   2.0)  Ci = complex<double>(AiryBiSeries(x),AiryAiSeries(x));
   else if (x <  10.0)  Ci = complex<double>(AiryBiSeries(x),AiryAiSteed(x));
   else                 Ci = complex<double>(AiryBiPosAsym(x),AiryAiPosAsym(x));

   return Ci;
}

inline complex<double> AiryCiPrime(double x)
{
/* *** Calculate the Airy function Ci'(z)
   *** July 26, 2002, written by CB */

   complex<double> CiP;
   if      (x < - 6.2)  CiP = AiryCiPrimeNegAsym(x);
   else if (x <   2.0)  CiP = complex<double>(AiryBiPrimeSeries(x),AiryAiPrimeSeries(x));
   else if (x <  10.0)  CiP = complex<double>(AiryBiPrimeSeries(x),AiryAiPrimeSteed(x));
   else                 CiP = complex<double>(AiryBiPrimePosAsym(x),AiryAiPrimePosAsym(x));

   return CiP;
}


/* *** FUNCTIONS: POSITION *** */

/* Switch Between Cartesian and Polar Coordinates */
	inline double PolarToX(double rh, double Ph)
	{
		return rh * cos(Ph);
	}
	inline double PolarToY(double rh, double Ph)
	{
		return rh * sin(Ph);
	}
	inline double XYToRho(double xx, double yy)
	{
		return sqrt(xx * xx + yy * yy);
	}

/* Calculate Dimensionless Position as a Function of Time and Emission Angles - Cartesian Coordinates */
	inline double xPos(double th, double ph, double tau)
	{
		return sin(th) * cos(ph - tau) * sin(tau);
	}
	inline double yPos(double th, double ph, double tau)
	{
		return sin(th) * sin(ph - tau) * sin(tau);
	}
	inline double zPos(double th, double tau)
	{
		return tau * cos(th);
	}

/* Calculate Dimensionless Position as a Function of Time and Emission Angles - Lateral Distance */
	inline double rhoPos(double th, double tau)
	{
		return sin(th) * fabs(sin(tau));
	}

/* *** FUNCTIONS: ANGLES (DIMENSIONLESS) *** */

/* Switch Between Cartesian and Polar Coordinates */
	inline double XYToPhi(double xx, double yy)
	{
		complex<double> Aux = complex<double>(xx, yy);
		return arg(Aux);
	}

/* Calculate Final Angle as a Function of Time and Emission Angle (and reverse) */
	inline double phiToPhi(double ph, double tau)
	{
		double Aux = ph - tau;
		return Aux - 2.0 * Pi * floor(Aux / (2.0 * Pi));
	}
	inline double PhiTophi(double Ph, double tau)
	{
		double Aux = Ph + tau;
		return Aux - 2.0 * Pi * floor(Aux / (2.0 * Pi));
	}

/* Calculate Emission Angle as a Function of Time of Flight and Final Position */
	inline double thetaAngle(double zz, double tau)
	{
		return acos(zz / tau);
	}
	inline double phiAngle(double xx, double yy, double tau)
	{
	/* Final Polar Angle */
		double Aux = XYToPhi(xx, yy);

	/* Correct For Sign of sin(tau) */
		if (sin(tau) < 0.0)  Aux += Pi;

	/* Add Contribution Along Trajectory and Bring Into Range */
		return PhiTophi(Aux, tau);
	}

/* *** FUNCTIONS: VELOCITY (dimensionless) *** */

/* Calculate Velocity as a Function of Time and Emission Angles */
	inline double xVel(double th, double ph, double tau)
	{
		return sin(th) * cos(ph - 2.0 * tau);
	}
	inline double yVel(double th, double ph, double tau)
	{
		return sin(th) * sin(ph - 2.0 * tau);
	}
	inline double zVel(double th, double tau)
	{
		return cos(th);
	}

/* Calculate Initial Velocity as a Function of Time of Flight and Final Position */
	inline double xVelInit(double xx, double yy, double tau)
	{
		return xx * cot(tau) - yy;
	}
	inline double yVelInit(double xx, double yy, double tau)
	{
		return yy * cot(tau) + xx;
	}
	inline double zVelInit(double zz, double tau)
	{
		return zz / tau;
	}

/* Calculate Final Velocity as a Function of Time of Flight and Final Position */
	inline double xVelFinal(double xx, double yy, double tau)
	{
		return xx * cot(tau) + yy;
	}
	inline double yVelFinal(double xx, double yy, double tau)
	{
		return yy * cot(tau) - xx;
	}
	inline double zVelFinal(double zz, double tau)
	{
		return zz / tau;
	}

/* *** FUNCTIONS: MOMENTUM (dimensionless) *** */

/* Calculate Canonical Momentum as a Function of Time and Emission Angles */
	inline double xMom(double th, double ph, double tau)
	{
		return sin(th) * cos(tau) * cos(ph - tau);
	}
	inline double yMom(double th, double ph, double tau)
	{
		return sin(th) * cos(tau) * sin(ph - tau);
	}
	inline double zMom(double th, double tau)
	{
		return cos(th);
	}

/* Calculate Initial Momentum as a Function of Time of Flight and Final Position */
	inline double xMomInit(double xx, double yy, double tau)
	{
		return xx * cot(tau) - yy;
	}
	inline double yMomInit(double xx, double yy, double tau)
	{
		return yy * cot(tau) + xx;
	}
	inline double zMomInit(double zz, double tau)
	{
		return zz / tau;
	}

/* Calculate Final Momentum as a Function of Time of Flight and Final Position */
	inline double xMomFinal(double xx, double yy, double tau)
	{
		return xx * cot(tau);
	}
	inline double yMomFinal(double xx, double yy, double tau)
	{
		return yy * cot(tau);
	}
	inline double zMomFinal(double zz, double tau)
	{
		return zz / tau;
	}

/* *** FUNCTIONS: ACTION, ENERGY, CROSS SECTION (dimensionless) *** */

/* Calculate Time-Dependent Action Integral For Final Position (real time) */
	inline double Action(double rsq, double zz, double eps, double tau)
	{
		return eps * (rsq * cot(tau) + zz * zz / tau);
	}

/* Calculate Time-Dependent Action Integral For Final Position (complex time) */
	inline complex<double> ComplexAction(double rsq, double zz, double eps, complex<double> tau)
	{
		return eps * (rsq * (cos(tau) / sin(tau)) + zz * zz / tau);
	}

/* Calculate Energy For Final Position, Time of Flight (real time) */
	inline double Energy(double rsq, double zz, double eps, double tau)
	{
		double av = 1.0 / sin(tau);
		double viz = zz / tau;
		return eps * (av * av * rsq + viz * viz);
	}

/* Calculate Energy For Final Position, Time of Flight (complex time) */
	inline complex<double> ComplexEnergy(double rsq, double zz, double eps, complex<double> tau)
	{
		complex<double> av = 1.0 / sin(tau);
		complex<double> viz = zz / tau;
		return eps * (av * av * rsq + viz * viz);
	}

/* Calculate Reduced Action Integral For Final Position (real action) */
	inline double ReducedAction(double rsq, double zz, double eps, double tau)
	{
	/* Check Numerics, If Desired */
		if (CheckNumerics == true)
		{
			double NumericalEnergy = Energy(rsq, zz, eps, tau);
			double Deviation = fabs(eps - NumericalEnergy);

			if (Deviation > (Accu * eps))
			{
				cerr << "\nWARNING:  Inaccurate Time-Energy Pair Encountered\n";
				cerr << "at (" << sqrt(rsq) << ", " << zz << ") for real time tau = " << tau << "\n";
				cerr << "Relative deviation of energy: " << Deviation / eps << "\n\n";
			}
		}

	/* Otherwise: Calculate Reduced Action From Time of Flight */
		return Action(rsq, zz, eps, tau) + eps * tau;
	}

/* Calculate Reduced Action Integral For Final Position (complex action) */
	inline complex<double> ComplexReducedAction(double rsq, double zz, double eps, complex<double> tau)
	{
	/* Check Numerics, If Desired */
		if (CheckNumerics == true)
		{
			complex<double> NumericalEnergy = ComplexEnergy(rsq, zz, eps, tau);
			double DeviationSquare = norm(eps - NumericalEnergy);

			if (DeviationSquare > Accu * Accu * eps * eps)
			{
				cerr << "\nWARNING:  Inaccurate Time-Energy Pair Encountered\n";
				cerr << "at (" << sqrt(rsq) << ", " << zz << ") for complex time tau = " << tau << "\n";
				cerr << "Relative deviation of energy: " << sqrt(DeviationSquare) / eps << "\n\n";
			}
		}

	/* Otherwise: Calculate Reduced Action From Time of Flight */
		return ComplexAction(rsq, zz, eps, tau) + eps * tau;
	}


/* Second Temporal Derivative of the Action Functional (real time) */
	inline double D2Action(double rsq, double zz, double eps, double tau)
	{
		double Aux = 1.0 / sin(tau);

		return 2.0 * eps * (Aux * Aux * Aux * rsq * cos(tau) + zz * zz / (tau * tau * tau));
	}

/* Second Temporal Derivative of the Action Functional (complex time) */
	inline complex<double> ComplexD2Action(double rsq, double zz, double eps, complex<double> tau)
	{
		complex<double> Aux = 1.0 / sin(tau);

		return 2.0 * eps * (Aux * Aux * Aux * rsq * cos(tau) + zz * zz / (tau * tau * tau));
	}

/* Third Temporal Derivative of the Action Functional */
	inline double D3Action(double rsq, double zz, double eps, double tau)
	{
		double sine = sin(tau);
		double zsq = zz * zz;
		double tausq = tau * tau;
		double sinsq = sine * sine;
		double secsq = 1.0 / sinsq;

		return -2.0 * eps * (secsq * secsq * rsq * (3.0 - 2.0 * sinsq) + 3.0 * zsq / (tausq * tausq));
	}

/* Classical Density (real time) */
/* (in free-particle units of (m k / 4 pi epsilon hbar^2)^2 */
	inline double Density(double rsq, double zz, double eps, double tau)
	{
		double sine = sin(tau);
		double zsq = zz * zz;
		double tausq = tau * tau;
		double sin3p = sine * sine * sine;

		return fabs((tausq * sine) / (rsq * tausq * tau * cos(tau) + zsq * sin3p));
	}

/* Classical Density (complex time) */
	inline complex<double> ComplexDensity(double rsq, double zz, double eps, complex<double> tau)
	{
		complex<double> sine = sin(tau);
		double zsq = zz * zz;
		complex<double> tausq = tau * tau;
		complex<double> sin3p = sine * sine * sine;

		return -tausq * sine / (rsq * tausq * tau * cos(tau) + zsq * sin3p);
	}

/* Classical Current Densities */
	inline double xCurrDens(double xx, double yy, double zz, double eps, double tau)
	{
		double rsq = xx * xx + yy * yy;
		return xVelFinal(xx, yy, tau) * Density(rsq, zz, eps, tau);
	}

	inline double yCurrDens(double xx, double yy, double zz, double eps, double tau)
	{
		double rsq = xx * xx + yy * yy;
		return yVelFinal(xx, yy, tau) * Density(rsq, zz, eps, tau);
	}

	inline double zCurrDens(double xx, double yy, double zz, double eps, double tau)
	{
		double rsq = xx * xx + yy * yy;
		return zVelFinal(zz, tau) * Density(rsq, zz, eps, tau);
	}

/* *** ROOT-FINDING ROUTINES (dimensionless) *** */

/* Find Left Boundary of Time Interval Containing Solutions */
	inline double TMinus(double zz, double eps)
	{
	/* Minimal Time (for motion along axis) */
		return fabs(zz);
	}

/* Find Time of Flight For Minimum Energy in Interval (k-1) * Pi < tau < k * Pi */
/*  (Note - fixed Oct 2011) */

	double MinimalTime(double rsq, double zz, double eps, int k)
	{
	/* Determine Time of Flight For Minimum Energy - Newton Method */
		int Counter = 0;

	/* Interval Limits */
		double tMax = k * Pi;
		double tMin = tMax - Pi;

	/* Starting Value:  Center of Interval */
		double t0 = tMax - PiOver2;

		double zsq = zz * zz;

	/* If z=0, Work is Done */
		if (zsq > 0)
		{
		/* Calculate Approximate Time of Flight */
			double DeltaT = atan(t0 * pow(rsq / zsq, 1.0/3.0));

			t0 = k * Pi - DeltaT;

		/* Newton Loop */
			while (fabs(DeltaT) > Accu)
			{
				DeltaT = - D2Action(rsq, zz, eps, t0) / D3Action(rsq, zz, eps, t0);
				++Counter;
				t0 += DeltaT;

			/* Check for Out of Interval - Convergence Failed? */
				if ((t0 > tMax) || (t0 < tMin))
				{
					cerr << "\nERROR:  Newton Iteration Did Not Converge In Function MinimalTime(...)\n\n";
					exit(1);
				}

			/* Failure To Converge in 30 Iterations - Numerical Accuracy Problem? */
				if (Counter > 30)
				{
					cerr << "\nWARNING:  Newton Iteration Failed Accuracy Goal In Function MinimalTime(...)\n";
					cerr <<   "          (relative accuracy: " << fabs(DeltaT) / Pi << ")\n\n";
					DeltaT = 0.0;
				}
			}
		}

	/* Return Minimum Time */
		return t0;
	}

/* Find Minimum of Energy in Interval (k-1)*Pi < tau < k*Pi */
	double MinimalEnergy(double rsq, double zz, double eps, int k)
	{
	/* Determine Time of Flight For Minimum Energy - Newton Method */
		double t0 = MinimalTime(rsq, zz, eps, k);

	/* Return Minimum Energy */
		return Energy(rsq, zz, eps, t0);
	}

/* Find Solutions for Time of Flight At Given Energy in Interval (k-1)*Pi < tau < k*Pi */
/* (the boolean variable fast determines whether the "fast" or "slow" solution is found) */
	double TimeOfFlight(double rsq, double zz, double eps, int k, bool fast)
	{
	/* Determine Time of Flight For Minimum Energy - Newton Method */
		double DeltaT = 1.0;
		int Counter = 0;

	/* Interval Limits */
		double tMax = k * Pi;
		double tMin = tMax - Pi;

	/* Select Appropriate Starting Value */
		double t0 = asin(sqrt(rsq));
		if (fast == true)
		{
			t0 += tMin;
		}
		else
		{
			t0 = tMax - t0;
		}

	/* Newton Loop */
		while (fabs(DeltaT) > Accu)
		{
		/* Determine Ratio (Energy - Epsilon) / D2Action */
			double sine = sin(t0);
			double cosine = cos(t0);
			double sin3p = sine * sine * sine;
			double t03p = t0 * t0 * t0;

			DeltaT = 0.5 * (rsq * t03p * sine + zz * zz * sin3p * t0 - sin3p * t03p) / (rsq * t03p * cosine + zz * zz * sin3p);

			++Counter;
			t0 += DeltaT;

		/* Check for Out of Interval - Convergence Failed? */
			if ((t0 > tMax) || (t0 < tMin))
			{
				cerr << "\nERROR:  Newton Iteration Did Not Converge In Function TimeOfFlight(...)\n\n";
				exit(1);
			}

		/* Failure To Converge in 30 Iterations - Numerical Accuracy Problem? */
			if (Counter > 30)
			{
				cerr << "\nWARNING:  Newton Iteration Failed Accuracy Goal In Function TimeOfFlight(...)\n";
				cerr <<   "          (relative accuracy: " << fabs(DeltaT) / Pi << ")\n\n";
				DeltaT = 0.0;
			}
		}

	/* Return Time of Flight */
		return t0;
	}

/* Find Complex Solutions for Time of Flight At Given Energy in Interval (k-1)*Pi < Re[tau] < k*Pi */
/* (the boolean variable upper determines whether solution in upper or lower complex plane is found) */
	complex<double> ComplexTimeOfFlight(double rsq, double zz, double eps, int k, bool upper)
	{
	/* Determine Time of Flight - Newton Method */
		complex<double> DeltaT = complex<double> (1.0, 0.0);
		int Counter = 0;

	/* Real Interval Limits */
		double tMax = k * Pi;
		double tMin = tMax - Pi;

	/* Select Appropriate Starting Value */
	/* Step 1: Find Real Time of Flight for Minimum of Energy */
		double tReal = MinimalTime(rsq, zz, eps, k);

	/* Step 2: Find Minimum of Energy Itself */
		double EMin = Energy(rsq, zz, eps, tReal);

	/* Step 3: Find Approximate Imaginary Part of TOF */
		double tImag = sqrt(2.0 * (eps - EMin) / D3Action(rsq, zz, eps, tReal));

	/* Step 4: Pick Proper Sign */
		if (upper == false)  tImag = -tImag;
		complex<double> t0 = complex<double> (tReal, tImag);

	/* Newton Loop */
		while (norm(DeltaT) > Accu * Accu)
		{
			DeltaT = (ComplexEnergy(rsq, zz, eps, t0) - eps) / ComplexD2Action(rsq, zz, eps, t0);
			++Counter;
			t0 += DeltaT;

		/* Check for Out of Interval - Convergence Failed? */
			if ((real(t0) < tMin) || (real(t0) > tMax))
			{
				cerr << "\nERROR:  Newton Iteration Did Not Converge In Function ComplexTimeOfFlight(...)\n\n";
				exit(1);
			}

		/* Failure To Converge in 30 Iterations - Numerical Accuracy Problem? */
			if (Counter > 30)
			{
				cerr << "\nWARNING:  Newton Iteration Failed Accuracy Goal In Function ComplexTimeOfFlight(...)\n";
				cerr <<   "          (relative accuracy: " << sqrt(norm(DeltaT / Pi)) << ")\n\n";
				DeltaT = 0.0;
			}
		}

	/* Return Time of Flight */
		return t0;
	}


/* *** MAIN PROGRAM ************************************************************** */

int _tmain(int argc, _TCHAR* argv[])
{
	/* Opening Message */
	cerr << "MAGNETIC CURRENT PROFILE.c++ - Semiclassical Current in a Homogeneous Magnetic Field\n"
		 << "(Version " << VersionSpec << ", written " << VersionDate << ")\n\n";

/* *** STEP 1: CALCULATE A DENSITY PROFILE ***/

/* Variables */
	long j, jj;
	double x, y, Aux, Aux2;
	complex<double> Sum, SumP, CAux;

/* Start Message */
	cerr << "Calculating Density Profile: \n";

/* For Analysis Purposes */
	double DensitySum = 0;
	double DensityMax = 0;
	double DensityMin = 0;
	double CurrentSum = 0;

/* Reserve Field for Density Information */
	double *DensityClassical = new double[xSamples];
	double *DensitySC = new double[xSamples];
	double *DensityUniform = new double[xSamples];
	double *DensityQuantum = new double[xSamples];
	double *CurrentDensityClassical = new double[xSamples];
	double *CurrentDensitySC = new double[xSamples];
	double *CurrentDensityUniform = new double[xSamples];
	double *CurrentDensityQuantum = new double[xSamples];

/* Just to avoid trouble with rho = 0 */
	y = 1e-9;

/* Establish Indices of Interval with Possible Solutions */
	int kMin = (int)ceil(TMinus(zDisp, Epsilon) / Pi);
	int kMax = (int)MaxOrbit;

/* Reserve Fields For Classical Flags, Solutions */
	bool *Exists = new bool[kMax - kMin + 1];
	double *TOF  = new double[2 * (kMax - kMin + 1)];
	complex<double> *ComplexTOF  = new complex<double>[2 * (kMax - kMin + 1)];

/* Reserve Array For 1D Green Functions */
	complex<double> *Green1DArray = new complex<double>[MaxLevel];
	complex<double> *Green1DPrimeArray = new complex<double>[MaxLevel];

/* Prepare Array of 1D Free Particle Green Function */

/* Landau Level Counter */
	long Level = 0;

/* Find Number of Open Channels */
	long Channel = (long)ceil((Epsilon - 1.0) / 2.0);

/* Green Function - Open Channels */
	while (Level < Channel)
	{
	/* Longitudinal Energy Parameter */
		Aux = sqrt(Epsilon * (Epsilon - (2 * Level + 1)));

	/* One-Dimensional Green Function in Field Direction, and its Derivative */
		Green1DPrimeArray[Level] = polar(1.0, 2.0 * Aux * fabs(zDisp));
		Green1DArray[Level] = -ii * Green1DPrimeArray[Level] / Aux;

	/* Apply Sign Correction, If Necessary */
		if (zDisp < 0) Green1DPrimeArray[Level] *= (-1);

	/* Increment Landau Level */
		++Level;
	}

/* Green Function - Closed Channels */
	do
	{
	/* Longitudinal Energy Parameter */
		Aux = sqrt(Epsilon * ((2 * Level + 1) - Epsilon));

	/* One-Dimensional Green Function in Field Direction, Derivative */
		Aux2 = exp(-2.0 * Aux * fabs(zDisp));
		Green1DArray[Level] = - Aux2 / Aux;
		Green1DPrimeArray[Level] = Aux2;

	/* Apply Sign Correction, If Necessary */
		if (zDisp < 0) Green1DPrimeArray[Level] *= (-1);

	/* Increment Landau Level */
		++Level;
	}
	while ((Aux2 > Accu) && (Level < MaxLevel));

/* Indicate Progress and Convergence Failure */
	if (Level == MaxLevel) 
		cerr << "WARNING:  Convergence not achieved in quantum series solution \n";

/* Calculate Densities and Currents - Loop Through Radial Distance */
	for (jj = 0; jj < xSamples; ++jj)
	{
		x = xMin + jj * (xMax - xMin) / (xSamples - 1);

	/* Radial Distance Parameters */
		double rSquare = x * x + y * y;
		double r = sqrt(rSquare);

	/* TRAJECTORY-BASED METHODS */

	/* Establish Existence of Real Solutions, Times of Flight, Average Action */
		int NumberOfSolutions = 0;

		for (j = kMin; j <= kMax; ++j)
		{
			if (MinimalEnergy(rSquare, zDisp, Epsilon, j) < Epsilon)
			{
			/* *** REAL (CLASSICAL) TRAJECTORIES *** */

				Exists[j - kMin] = true;
				NumberOfSolutions += 2;

			/* Time of Flight - Fast Track */
				TOF[2 * (j - kMin)]      = TimeOfFlight(rSquare, zDisp, Epsilon, j, true);

			/* Time of Flight - Slow Track */
				TOF[2 * (j - kMin) + 1]  = TimeOfFlight(rSquare, zDisp, Epsilon, j, false);
			}
			else
			{
			/* *** COMPLEX (TUNNELING) TRAJECTORIES *** */
				Exists[j - kMin] = false;

			/* Time of Flight - Solution in Upper Half Plane */
				ComplexTOF[2 * (j - kMin)] = ComplexTimeOfFlight(rSquare, zDisp, Epsilon, j, true);

			/* Time of Flight - Solution in Lower Half Plane */
				ComplexTOF[2 * (j - kMin) + 1] = conj(ComplexTOF[2 * (j - kMin)]);
			}
		}

	/* Classical Density */
		double ClassicalDensity = 0.0;
		double ClassicalCurrentDensity = 0.0;

	/* Semiclassical Density */
		int Maslov;
		double Phase;
		complex<double> WaveFunction = complex<double>(0.0, 0.0);
		complex<double> WaveFunctionPrime = complex<double>(0.0, 0.0);

	/* Uniform Approximation */
		double PhaseFast, PhaseSlow;
		double AmplitudeFast, AmplitudeSlow;
		double AiryAmplitude, AiryPrimeAmplitude, AiryArg;
		complex<double> Amplitude;
		complex<double> UniformWaveFunction = complex<double>(0.0, 0.0);
		complex<double> UniformWaveFunctionPrime = complex<double>(0.0, 0.0);

		for (j = kMin; j <= kMax; ++j)
		{
			if (Exists[j - kMin] == true)
			{
			/* *** CLASSICALLY ALLOWED MOTION *** */
			/* *** FAST TRACK *** */

			/* Time Of Flight */
				double tFast = TOF[2 * (j - kMin)];

			/* Maslov Index */
				Maslov = 2 * (j - 1);

			/* Phase of the Wave Function */
				PhaseFast = ReducedAction(rSquare, zDisp, Epsilon, tFast);
				Phase = PhaseFast - 0.5 * Pi * Maslov;

			/* Amplitude of the Wave Function */
				double rho = Density(rSquare, zDisp, Epsilon, tFast);
				AmplitudeFast = sqrt(rho);
			
			/* Classical Particle Velocity in z Direction */
				double zVelFast = zDisp / tFast;
			
			/* Classical Density */
				ClassicalDensity += rho;
				ClassicalCurrentDensity += zVelFast * rho;

			/* Add To Wave Function */
				CAux = polar(AmplitudeFast, Phase);
				WaveFunction -= CAux;
				WaveFunctionPrime -= zVelFast * CAux;

			/* *** SLOW TRACK *** */

			/* Time Of Flight */
				double tSlow = TOF[2 * (j - kMin) + 1];

			/* Maslov Index */
				Maslov = 2 * j - 1;

			/* Phase of the Wave Function */
				PhaseSlow = ReducedAction(rSquare, zDisp, Epsilon, tSlow);
				Phase = PhaseSlow - 0.5 * Pi * Maslov;

			/* Amplitude of the Wave Function */
				rho = Density(rSquare, zDisp, Epsilon, tSlow);
				AmplitudeSlow = sqrt(rho);

			/* Classical Particle Velocity in z Direction */
				double zVelSlow = zDisp / tSlow;

			/* Classical Density */
				ClassicalDensity += rho;
				ClassicalCurrentDensity += zVelSlow * rho;

			/* Add To Wave Function */
				CAux = polar(AmplitudeSlow, Phase);
				WaveFunction -= CAux;
				WaveFunctionPrime -= zVelSlow * CAux;

			/* Uniform Approximation */
		
			/* Calculate Airy Functions */
				Phase = ((double)j - 0.25) * Pi + 0.5 * (PhaseFast + PhaseSlow);
				AiryArg = - pow(0.75 * (PhaseSlow - PhaseFast), 2.0 / 3.0);
				double Ai = AiryAi(AiryArg);
				double AiPrime = AiryAiPrime(AiryArg);
				double AirySqrt = sqrt(Pi * sqrt(-AiryArg));

			/* Uniform Approximation for Wave Function */

			/* Find Weights of Airy Functions */
				AiryAmplitude = AirySqrt * (AmplitudeFast + AmplitudeSlow);
				AiryPrimeAmplitude = (Pi / AirySqrt) * (AmplitudeFast - AmplitudeSlow);

			/* Add To Uniform Wave Function */
				UniformWaveFunction += polar(AiryAmplitude, Phase) * Ai;
				UniformWaveFunction += polar(AiryPrimeAmplitude, Phase - PiOver2) * AiPrime;

			/* Uniform Approximation for z-Derivative of Wave Function */

			/* Adjust Amplitudes */
				AmplitudeSlow *= zVelSlow;
				AmplitudeFast *= zVelFast;

			/* Find Weights of Airy Functions */
				AiryAmplitude = AirySqrt * (AmplitudeFast + AmplitudeSlow);
				AiryPrimeAmplitude = (Pi / AirySqrt) * (AmplitudeFast - AmplitudeSlow);

			/* Add To Uniform Wave Function */
				UniformWaveFunctionPrime += polar(AiryAmplitude, Phase + PiOver2) * Ai;
				UniformWaveFunctionPrime += polar(AiryPrimeAmplitude, Phase) * AiPrime;
			}
		else
			{
			/* *** TUNNELING TRAJECTORIES *** */

			/* *** LOWER HALF OF COMPLEX PLANE *** */

			/* Time Of Flight */
				complex<double> ComplexT = ComplexTOF[2 * (j - kMin) + 1];

			/* Phase of the Wave Function */
				CAux = ComplexReducedAction(rSquare, zDisp, Epsilon, ComplexT);
				Phase = (j + 0.5) * Pi + real(CAux);
				double Decay = imag(CAux);

			/* Amplitude of the Wave Function */
				Amplitude = sqrt(ComplexDensity(rSquare, zDisp, Epsilon, ComplexT));

			/* Complex Velocity in z Direction */
				complex<double> zVelComplex = zDisp / ComplexT;

			/* Add To Semiclassical Wave Function, If Desired */
				if (IncludeTunneling == true)
				{
					CAux = Amplitude * polar(exp(-Decay), Phase);
					WaveFunction -= CAux;
					WaveFunctionPrime -= zVelComplex * CAux;
				}

			/* Uniform Approximation */

			/* Calculate Airy Functions */
				AiryArg = pow(1.5 * Decay, 2.0 / 3.0);
				double Ai = AiryAi(AiryArg);
				double AiPrime = AiryAiPrime(AiryArg);
				double AirySqrt = sqrt(2.0 * Pi * sqrt(AiryArg));

			/* Uniform Approximation for Wave Function */

			/* Find Weights of Airy Functions */
				AiryAmplitude = AirySqrt * (real(Amplitude) + imag(Amplitude));
				AiryPrimeAmplitude = (2.0 * Pi / AirySqrt) * (imag(Amplitude) - real(Amplitude));

			/* Add To Uniform Wave Function */
				Phase -= 0.75 * Pi;
				UniformWaveFunction += polar(AiryAmplitude, Phase) * Ai;
				UniformWaveFunction += polar(AiryPrimeAmplitude, Phase - PiOver2) * AiPrime;

			/* Uniform Approximation for z-Derivative of Wave Function */

			/* Adjust Amplitude */
				Amplitude *= zVelComplex;

			/* Find Weights of Airy Functions */
				AiryAmplitude = AirySqrt * (real(Amplitude) + imag(Amplitude));
				AiryPrimeAmplitude = (2.0 * Pi / AirySqrt) * (imag(Amplitude) - real(Amplitude));

			/* Add To Uniform Wave Function */
				UniformWaveFunctionPrime += polar(AiryAmplitude, Phase + PiOver2) * Ai;
				UniformWaveFunctionPrime += polar(AiryPrimeAmplitude, Phase) * AiPrime;
			}
		}

	/* EXACT QUANTUM SOLUTION */

	/* Argument of Laguerre Polynomial */
		double Radial = 2.0 * Epsilon * rSquare;

	/* Initial Recursion Values For Laguerre Polynomials */
		double l0  = 1.0;
		double l1  = 1.0 - Radial;
		double l2;

	/* Calculate Green Function via Series */
		Sum = 0.0;
		SumP = 0.0;

		for (long n = 0; n < Level; ++n)
		{
		/* Add Term To Sum */
			Sum += l0 * Green1DArray[n];
			SumP += l0 * Green1DPrimeArray[n];
	
		/* Recursion Step For Laguerre Polynomials */			
			l2  = ((2 * n + 3 - Radial) * l1 - (n + 1) * l0) / (n + 2);
			l0  = l1;
			l1  = l2;
		}

	/* Apply Prefactor */
		Aux = exp(-Radial);
			
	/* Store Density, Current Density Vector */
		double QuantumDensity          =  4.0 * Aux * Epsilon * Epsilon * norm(Sum);
		double QuantumCurrentDensity   =  4.0 * Aux * Epsilon * imag(conj(Sum) * SumP);

	/* Integrated Densities and Currents, If Desired */
		if (IntegratedDensity == true)
		{
			Aux = 2.0 * Pi * fabs(x);
		}
		else
		{
			Aux = 1.0;
		}

	/* Calculate Desired Result(s) */
		DensityClassical[jj] = Aux * ClassicalDensity;
		DensitySC[jj] = Aux * norm(WaveFunction);
		DensityUniform[jj] = Aux * norm(UniformWaveFunction);
		DensityQuantum[jj] = Aux * QuantumDensity;

		CurrentDensityClassical[jj] = Aux * ClassicalCurrentDensity;
		CurrentDensitySC[jj] = Aux * real(conj(WaveFunction) * WaveFunctionPrime);
		CurrentDensityUniform[jj] = Aux * imag(conj(UniformWaveFunction) * UniformWaveFunctionPrime);
		CurrentDensityQuantum[jj] = Aux * QuantumCurrentDensity;

	/* Show Progress */
		cerr << ".";
	}

/* Clean-Up */
	delete[] Exists;
	delete[] TOF;
	delete[] ComplexTOF;

/* Notify User */
	cerr << " done.\n\n";

/* *** STEP 3: WRITE DATA FILE *** */

/* Open Image File in Binary Mode */
	ofstream OutputDensity (OutputFile.c_str(), (ios::out));

/* Check For Success */
	if (!OutputDensity)
	{
		cerr << "ERROR:  Could Not Create Density File. \n\n";
		exit(-1);
	}

	cerr << "Writing Data File ... ";

/* Format Numbers */
	OutputDensity.precision(6);

/* Write Brief Info Header */
	OutputDensity << "# SOURCE IN MAGNETIC FIELD - DENSITY AND CURRENT PROFILES \n";
	OutputDensity << "# (using dimensionless units)\n";
	OutputDensity << "# \n";
	OutputDensity << "# Energy parameter: " << Epsilon << "\n";
	OutputDensity << "# Displacement: " << zDisp << "\n";
	OutputDensity << "# Integrated Current: " << IntegratedDensity << "\n";
	OutputDensity << "# Includes Tunneling: " << IncludeTunneling << "\n";
	OutputDensity << "# \n";

	OutputDensity << "# Output Structure: \n";
	OutputDensity << "# [1] Radius rho \n";
	OutputDensity << "# [2] Classical Density \n";
	OutputDensity << "# [3] Semiclassical Density \n";
	OutputDensity << "# [4] Quantum Density \n";
	OutputDensity << "# [5] Uniform Density \n";
	OutputDensity << "# [6] Classical Current (along field) \n";
	OutputDensity << "# [7] Semiclassical Current \n";
	OutputDensity << "# [8] Uniform Current \n";
	OutputDensity << "# [9] Quantum Current \n";
	OutputDensity << "# \n";

/* Write Position, Current Info */
	for (jj = 0; jj < xSamples; ++jj)
	{
		x = xMin + jj * (xMax - xMin) / (xSamples - 1);

		OutputDensity << fixed 
			<< x << " \t"
			<< DensityClassical[jj] << " \t"
			<< DensitySC[jj] << " \t"
			<< DensityUniform[jj] << " \t"
			<< DensityQuantum[jj] << " \t"
			<< CurrentDensityClassical[jj] << " \t"
			<< CurrentDensitySC[jj] << " \t"
			<< CurrentDensityUniform[jj] << " \t"
			<< CurrentDensityQuantum[jj] << "\n";
	}
	
/* Clean Up Fields */
	delete[] DensityClassical;
	delete[] DensitySC;
	delete[] DensityUniform;
	delete[] DensityQuantum;
	delete[] CurrentDensityClassical;
	delete[] CurrentDensitySC;
	delete[] CurrentDensityUniform;
	delete[] CurrentDensityQuantum;

	delete[] Green1DArray;
	delete[] Green1DPrimeArray;

/* Write File And Close */
	OutputDensity.close();

	cerr << " done.\n\n";

/* *** CLEAN-UP *** */
	return 0;
}


