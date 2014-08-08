/*****************************************************************************************

  Semiclassical Electron Density in a Purely Magnetic Field
  (using dimensionless parameters - free particle units)

  Written by:  Christian Bracher
  Version   :  1.02
  Date      :  December 7, 2011

*****************************************************************************************/

/* *** INITIALIZATION ROUTINES *** */
#include "stdafx.h"
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string>
#include<complex>

using namespace std;

/* Version Info */
	const char *VersionSpec = "1.02";
	const char *VersionDate = "December 7, 2011";

/* Bitmap Template */
	#pragma pack(push)
	#pragma pack(1)
	struct BMPHeader
	{
		short int BMPSymbol;
		long FileLength;
		long Unknown1;
		long HeaderLength;
		long Unknown2;
		long Width;
		long Height;
		short int Bitplanes;
		short int BitsPerPixel;
		long Unknown3;
		long BitmapSize;
		long ResolutionX;
		long ResolutionY;
		long Unknown4;
		long Unknown5;
	};
	#pragma pack(pop)

/* *** GLOBAL VARIABLES *** */

/* Mathematical Constants */
	double Pi = 4.0 * atan(1.0);
	double PiOver2 = Pi / 2.0;
	complex<double> ii = complex<double>(0.0,1.0);

/* *** Physical Parameters *** */

/* Initial Electron Energy [in units of hbar * omega_L] */
	double Epsilon = 51.01;

/* *** Math Control Parameters *** */

/* Accuracy Goal */
	double Accu = 1.0e-10;

/* *** Output Parameters *** */

/* Observation Window [in m] */
	double xMin = 0.001;
	double xMax = 1.1;

	double zMin = -1.1;
	double zMax = 3.3;

/* Sampling Number */
  	int xSamples = 1000;
	int zSamples = 4000;

/* The output file name */
  	string BitmapFile = "At 51 (semiclassical density).bmp";

/* *** CLASSICAL CUTOFF *** */

/* Maximum Number of Cyclotron Orbits */
	long MaxOrbit = 250;


/* *** PROGRAM CONTROL *** */

/* NOTE:  Program will return image for first "true" pick! */

/* Return Classical Density? */
	bool ReturnClassical = false;

/* Return Semiclassical Approximation? */
	bool ReturnSemiclassical = true;

/* Include Tunneling Trajectories in Semiclassical Approximation? */
	bool IncludeTunneling = true;

/* Return Uniform Approximation? */
	bool ReturnUniform = true;

/* Integrated Density? */
	bool IntegratedDensity = true;

/* Check for Numerical Inaccuracies? */
	bool CheckNumerics = false;

/* *** Graphics Control Parameters *** */

/* Frame Parameters */
	int FrameWidth = 2;
	int TickWidth = 2;
	int TickLength = 10;

	unsigned char FrameRed = (unsigned char)0;
	unsigned char FrameGreen = (unsigned char)0;
	unsigned char FrameBlue = (unsigned char)0;

/* Tick Levels */
	double xTick = 0.5;
	double zTick = 0.5;

/* Density Parameters */
	double MaxDensity = 100.0 * Epsilon;

/* Color Saturation Control */
	double Gamma = 0.5;
	double SaturationOffset = 1.0;


/* *** MATHEMATICAL ROUTINES *** */

/* *** Trigonometric Functions */
	inline double cot(double x);

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

/* *** The Laguerre Polynomials L[n,m](x) */
	double Laguerre(int n, int m, double x);

inline double cot(double x)
{
/* *** Calculate cot(x)
   *** October 2011, written by CB */

/* Return Function Value */
	return tan(PiOver2 - x);
}

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
/* (using "free particle units" (m k / 2 pi hbar^2)^2 ) */
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

/* *** GRAPHICS PROCEDURES ************************************************************ */

/* Clear the Image Bitmap Template (White Background) */
	void ClearBMPImage(char BMPImage[], long BMPWidth, long BMPHeight)
	{
		long BitmapSize = 3 * BMPHeight * BMPWidth;
		for (int i = 0; i < BitmapSize; ++i)  BMPImage[i] = char(0xFF);
	}

/* Paint a Pixel At Absolute Position (x,y) Within The Bitmap */
	void PSet(char BMPImage[], long BMPWidth, long BMPHeight, long XPixel, long YPixel,
		unsigned char Red, unsigned char Green, unsigned char Blue)
	{
	/* Check For Canvas Boundaries */
		if ((XPixel < 0) || (YPixel < 0) ||
			(XPixel >= BMPWidth) || (YPixel >= BMPHeight))
		{
			cerr << "WARNING:  Bitmap Dimensions Exceeded. \n\n";
		}
		else
		{
		/* 'Color' the Specified Pixel */
			long PixelPos = 3 * (BMPWidth * YPixel + XPixel);

			BMPImage[PixelPos]     = Blue;
			BMPImage[PixelPos + 1] = Green;
			BMPImage[PixelPos + 2] = Red;
		}
	}

/* Determine Pixel Color from Numerical Value [0 ... 1] */
	void AssignColor(unsigned char *Red, unsigned char *Green, unsigned char *Blue, double Value)
	{
		if ((Value < 0) || (Value > 1))
		{
			cerr << "ERROR:  Color Assignment Failure. \n\n";
			exit(1);
		}

      *Red   = (unsigned char)floor(255.99 * (1.0 - sqrt(Value)));
      *Green = (unsigned char)floor(255.99 * (1.0 - Value));
      *Blue  = (unsigned char)floor(255.99 * (1.0 - Value * Value));
	}


/* *** MAIN PROGRAM ************************************************************** */

int _tmain(int argc, _TCHAR* argv[])
{
	/* Opening Message */
	cerr << "SEMICLASSICAL MAP.c++ - Semiclassical Energy Green Function in a Homogeneous Magnetic Field\n"
		 << "(Version " << VersionSpec << ", written " << VersionDate << ")\n\n";

/* *** STEP 1: CALCULATE A DENSITY PROFILE ***/

/* Variables */
	long i, j, jj;
	double x, y, z, Aux;
	complex<double> Sum;

/* Start Message */
	cerr << "Calculating Density Profile: ";

/* Reserve Field for Density Information */
	double *ParticleDensity = new double[xSamples * zSamples];

/* Just To Avoid Trouble With r = 0 */
	y = 1e-9;

/* Loop Through z Values */
	for (i = 0; i < zSamples; ++i)
	{
	/* Linear z Values */
		z = zMin + i * (zMax - zMin) / (zSamples - 1);

	/* Establish Indices of Interval with Possible Solutions */
		int kMin = (int)ceil(TMinus(z, Epsilon) / Pi);
		int kMax = (int)MaxOrbit;

	/* Reserve Fields For Flags, Solutions */
		bool *Exists = new bool[kMax - kMin + 1];
		double *TOF  = new double[2 * (kMax - kMin + 1)];
		complex<double> *ComplexTOF  = new complex<double>[2 * (kMax - kMin + 1)];

	/* Calculate Green Function As A Sum */
		for (jj = 0; jj < xSamples; ++jj)
		{
			x = xMin + jj * (xMax - xMin) / (xSamples - 1);

		/* Radial Distance Parameters */

			double rSquare = x * x + y * y;
			double r = sqrt(rSquare);

		/* Establish Existence of Real Solutions, Times of Flight, Average Action */
			int NumberOfSolutions = 0;

			for (j = kMin; j <= kMax; ++j)
			{
				if (MinimalEnergy(rSquare, z, Epsilon, j) < Epsilon)
				{
				/* *** REAL (CLASSICAL) TRAJECTORIES *** */

					Exists[j - kMin] = true;
					NumberOfSolutions += 2;

				/* Time of Flight - Fast Track */
					TOF[2 * (j - kMin)]      = TimeOfFlight(rSquare, z, Epsilon, j, true);

				/* Time of Flight - Slow Track */
					TOF[2 * (j - kMin) + 1]  = TimeOfFlight(rSquare, z, Epsilon, j, false);
				}
				else
				{
				/* *** COMPLEX (TUNNELING) TRAJECTORIES *** */
					Exists[j - kMin] = false;

				/* Time of Flight - Solution in Upper Half Plane */
					ComplexTOF[2 * (j - kMin)] = ComplexTimeOfFlight(rSquare, z, Epsilon, j, true);

				/* Time of Flight - Solution in Lower Half Plane */
					ComplexTOF[2 * (j - kMin) + 1] = conj(ComplexTOF[2 * (j - kMin)]);
				}
			}


		/* Classical Density */
			double ClassicalDensity = 0.0;

		/* Semiclassical Density */
			int Maslov;
			double Phase;
			complex<double> WaveFunction = complex<double>(0.0, 0.0);

		/* Uniform Approximation */
			double PhaseFast, PhaseSlow;
			double AmplitudeFast, AmplitudeSlow;
			double AiryAmplitude, AiryPrimeAmplitude, AiryArg;
			complex<double> Amplitude;
			complex<double> UniformWaveFunction = complex<double>(0.0, 0.0);


			for (j = kMin; j <= kMax; ++j)
			{
				if (Exists[j - kMin] == true)
				{
				/* *** CLASSICALLY ALLOWED MOTION *** */
				/* *** FAST TRACK *** */

				/* Time Of Flight */
					double t = TOF[2 * (j - kMin)];

				/* Maslov Index */
					Maslov = 2 * (j - 1);

				/* Phase of the Wave Function */
					PhaseFast = ReducedAction(rSquare, z, Epsilon, t);
					Phase = PhaseFast - 0.5 * Pi * Maslov;

				/* Amplitude of the Wave Function */
					double rho = Density(rSquare, z, Epsilon, t);
					AmplitudeFast = sqrt(rho);

				/* Classical Density */
					ClassicalDensity += rho;

				/* Add To Wave Function */
					WaveFunction -= polar(AmplitudeFast, Phase);

				/* *** SLOW TRACK *** */

				/* Time Of Flight */
					t = TOF[2 * (j - kMin) + 1];

				/* Maslov Index */
					Maslov = 2 * j - 1;

				/* Phase of the Wave Function */
					PhaseSlow = ReducedAction(rSquare, z, Epsilon, t);
					Phase = PhaseSlow - 0.5 * Pi * Maslov;

				/* Amplitude of the Wave Function */
					rho = Density(rSquare, z, Epsilon, t);
					AmplitudeSlow = sqrt(rho);

				/* Classical Density */
					ClassicalDensity += rho;

				/* Add To Wave Function */
					WaveFunction -= polar(AmplitudeSlow, Phase);

				/* Uniform Approximation (if required) */

					if (ReturnUniform == true)
					{
						double TwoThird = 2.0 / 3.0;
						double PhaseAvg = 0.5 * (PhaseFast + PhaseSlow);
						double PhaseDiff = PhaseSlow - PhaseFast;

					/* Arguments and Prefactors */
						Phase = ((double)j - 0.25) * Pi + PhaseAvg;
						AiryArg = - pow(0.75 * PhaseDiff, TwoThird);

						double ArgRoot = sqrt(-AiryArg);
						double AmplitudeSum = AmplitudeFast + AmplitudeSlow;
						double AmplitudeDiff = AmplitudeFast - AmplitudeSlow;
						double RootRoot = sqrt(Pi * ArgRoot);

						AiryAmplitude = RootRoot * AmplitudeSum;
						AiryPrimeAmplitude = Pi * AmplitudeDiff / RootRoot;

					/* Add To Uniform Wave Function */
						UniformWaveFunction += polar(AiryAmplitude, Phase) * AiryAi(AiryArg);
						UniformWaveFunction += polar(AiryPrimeAmplitude, Phase - PiOver2) * AiryAiPrime(AiryArg);
					}
				}
				else
				{
				/* *** TUNNELING TRAJECTORIES *** */

				/* *** LOWER HALF OF COMPLEX PLANE *** */

				/* Time Of Flight */
					complex<double> ComplexT = ComplexTOF[2 * (j - kMin) + 1];

				/* Phase of the Wave Function */
					Phase = (2 * j + 1) * PiOver2 + real(ComplexReducedAction(rSquare, z, Epsilon, ComplexT));
					double Decay = imag(ComplexReducedAction(rSquare, z, Epsilon, ComplexT));

				/* Amplitude of the Wave Function */
					Amplitude = sqrt(ComplexDensity(rSquare, z, Epsilon, ComplexT));

				/* Add To Semiclassical Wave Function, If Desired */
					if (IncludeTunneling == true)
					{
						WaveFunction -= Amplitude * polar(exp(-Decay), Phase);
					}

				/* Uniform Approximation (if required) */

					if (ReturnUniform == true)
					{
					/* Arguments and Prefactors */
						double TwoThird = 2.0 / 3.0;
						double ReAmp = real(Amplitude);
						double ImAmp = imag(Amplitude);

						AiryArg = pow(1.5 * Decay, TwoThird);
						double AmpSum = ReAmp + ImAmp;
						double AmpDiff = ImAmp - ReAmp;

						double ArgRoot = sqrt(AiryArg);
						double RootRoot = sqrt(2.0 * Pi * ArgRoot);

						AiryAmplitude = RootRoot * AmpSum;
						AiryPrimeAmplitude = 2.0 * Pi * AmpDiff / RootRoot;

					/* Add To Uniform Wave Function */
						Phase -= 0.75 * Pi;
						UniformWaveFunction += polar(AiryAmplitude, Phase) * AiryAi(AiryArg);
						UniformWaveFunction += polar(AiryPrimeAmplitude, Phase - PiOver2) * AiryAiPrime(AiryArg);
					}
				}
			}

		/* Return Integrated Radial Density, If Desired */
			if (IntegratedDensity == true)
			{
				Aux = 2.0 * Pi * r;
			}
			else
			{
				Aux = 1.0;
			}

		/* Calculate Desired Result(s) */
			if (ReturnClassical == true)
			{
				Aux *= ClassicalDensity;
			}
			else
			{
				if (ReturnSemiclassical == true)
				{
					Aux *= norm(WaveFunction);
				}
				else
				{
					if (ReturnUniform == true)
					{
						Aux *= norm(UniformWaveFunction);
					}
				}
			}

		/* Store Density of Wave Function */
			ParticleDensity[i * xSamples + jj] = Aux;
		}

	/* Show Progress */
		cerr << ".";

	/* Clean-Up */
		delete[] Exists;
		delete[] TOF;
		delete[] ComplexTOF;
	}

/* Notify User */
	cerr << " done.\n\n";

/* *** STEP 2: CREATE BITMAP FILE *** */

/* Open Image File in Binary Mode */
	ofstream Bitmap (BitmapFile.c_str(), (ios::out | ios::binary));

/* Check For Success */
	if (!Bitmap)
	{
		cerr << "ERROR:  Could Not Create BMP Image File. \n\n";
		exit(-1);
	}

/* Write Bitmap File Header */
	long   BMPWidth          = xSamples;
	long   BMPHeight         = zSamples;
	int    PrinterResolution = 600;
	double InchPerMeter      = 39.3701;
	long   Resolution        = (long)(PrinterResolution * InchPerMeter);
	long   BitmapLength      = 3 * BMPWidth * BMPHeight;

	BMPHeader Header;
	Header.BMPSymbol  = 0x4d42;
	Header.FileLength = BitmapLength + 0x36;
	Header.Unknown1 = 0x0;
	Header.HeaderLength = 0x36;
	Header.Unknown2 = 0x28;
	Header.Width = BMPWidth;
	Header.Height = BMPHeight;
	Header.Bitplanes = 0x1;
	Header.BitsPerPixel = 0x18;
	Header.Unknown3 = 0x0;
	Header.BitmapSize = BitmapLength;
	Header.ResolutionX = Resolution;
	Header.ResolutionY = Resolution;
	Header.Unknown4 = 0x0;
	Header.Unknown5 = 0x0;

	Bitmap.write((char *)&Header, 54);

/* Create Canvas */
	char *ImageTemplate = new char[BitmapLength];
	ClearBMPImage(ImageTemplate, BMPWidth, BMPHeight);

	unsigned char Red   = ' ';
	unsigned char Green = ' ';
	unsigned char Blue  = ' ';

/* *** STEP 3: PAINT CANVAS *** */

/* Saturation Value for Image */
	double Saturation;

/* User Message */
	cerr << "Writing Bitmap Graphics File (with gamma value " << Gamma << "): ";

/* Loop Over z Values */
	for (i = 0; i < zSamples; ++i)
	{
	/* Loop Over x Values */
		for (j = 0; j < xSamples; ++j)
		{
		/* Calculate Relative Density */
			Aux = ParticleDensity[i * xSamples + j] / MaxDensity;

		/* Apply Saturation Offset */
			Aux *= SaturationOffset;
			if (Aux > 1.0) Aux = 1.0;

		/* Sort Out Unreliable Results */
			if (_isnan(Aux) != 0)
			{
				Red = (unsigned char)255;
				Green = (unsigned char)0;
				Blue = (unsigned char)0;
			}
			else
			{
			/* Adjust Color Intensity */
				Saturation = pow(Aux, Gamma);

			/* Assign Coloring */
				AssignColor(&Red, &Green, &Blue, Saturation);
			}

		/* Paint Pixel (Note: Progress top-to-bottom) */
			PSet(ImageTemplate, BMPWidth, BMPHeight, j, zSamples - i - 1, Red, Green, Blue);
		}

	/* Progress Marker */
		cerr << ".";
	}

/* Clean Up Fields */
	delete[] ParticleDensity;

/* Create Image Frame */
	for (i = 0; i < FrameWidth; ++i)
	{
		for (j = 0; j < zSamples; ++j)
		{
			PSet(ImageTemplate, BMPWidth, BMPHeight,
				i , j, FrameRed, FrameGreen, FrameBlue);
			PSet(ImageTemplate, BMPWidth, BMPHeight,
				xSamples - i - 1 , j, FrameRed, FrameGreen, FrameBlue);
		}
		for (j = 0; j < xSamples; ++j)
		{
			PSet(ImageTemplate, BMPWidth, BMPHeight,
				j , i, FrameRed, FrameGreen, FrameBlue);
			PSet(ImageTemplate, BMPWidth, BMPHeight,
				j, zSamples - i - 1, FrameRed, FrameGreen, FrameBlue);
		}
	}

/* Create Ticks */
	int xMinTick = (int)ceil(xMin / xTick);
	int xMaxTick = (int)floor(xMax / xTick);

	int zMinTick = (int)ceil(zMin / zTick);
	int zMaxTick = (int)floor(zMax / zTick);

	for (int l = xMinTick; l <= xMaxTick; ++l)
	{
		double Center = (l * xTick - xMin) * (xSamples - 1) / (xMax - xMin);
		long Upper = long(Center + 0.5 * TickWidth + 0.5);
		if (Upper > xSamples - 1) Upper = xSamples - 1;
		long Lower = long(Center - 0.5 * TickWidth + 0.5);
		if (Lower < 0) Lower = 0;

		for (i = Lower; i <= Upper; ++i)
		{
			for (j=0; j < TickLength; ++j)
			{
				PSet(ImageTemplate, BMPWidth, BMPHeight,
					i , j, FrameRed, FrameGreen, FrameBlue);
			}
			for (j = zSamples - 1; j >= zSamples - TickLength; --j)
			{
				PSet(ImageTemplate, BMPWidth, BMPHeight,
					i , j, FrameRed, FrameGreen, FrameBlue);
			}
		}
	}

	for (int l = zMinTick; l <= zMaxTick; ++l)
	{
	/* Position in z Space */
		double zCenter = l * zTick;

	/* Position in w Space */
		double Center = zCenter;

	/* Position on Canvas */
	 	Center = (Center - zMin) * (zSamples - 1) / (zMax - zMin);

		long Upper = long(Center + 0.5 * TickWidth + 0.5);
		if (Upper > zSamples - 1) Upper = zSamples - 1;
		long Lower = long(Center - 0.5 * TickWidth + 0.5);
		if (Lower < 0) Lower = 0;

		for (i = Lower; i <= Upper; ++i)
		{
			for (j = 0; j < TickLength; ++j)
			{
				PSet(ImageTemplate, BMPWidth, BMPHeight,
					j , zSamples - i - 1, FrameRed, FrameGreen, FrameBlue);
			}
			for (j = xSamples - 1; j >= xSamples - TickLength; --j)
			{
				PSet(ImageTemplate, BMPWidth, BMPHeight,
					j , zSamples - i - 1, FrameRed, FrameGreen, FrameBlue);
			}
		}
	}

/* Write Bitmap And Close */
	Bitmap.write(ImageTemplate, BitmapLength);
	Bitmap.close();

	cerr << " done.\n\n";

/* *** CLEAN-UP *** */

	delete[] ImageTemplate;

	return 0;
}

