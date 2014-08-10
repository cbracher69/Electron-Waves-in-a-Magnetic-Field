/*****************************************************************************************

  Quantum Electron Density and Current in a Purely Magnetic Field
  Dimensionless Units

  Written by:  Christian Bracher
  Version   :  1.03 
  Date      :  November 20, 2011

  NOTE:  Build for MS Visual C++ 2010 Express

*****************************************************************************************/

/* *** INITIALIZATION ROUTINES *** */

#include "stdafx.h"

/* Required Libraries */
#include<fstream>
#include<iostream>
#include<string>
#include<complex>

using namespace std;

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

/* Density and Current Template */

	struct DensCurr
	{
		double Density;
		double CurrentField;
		double CurrentRadial;
		double CurrentAzimuthal;
	};

/* Version Info */
	const char *VersionSpec = "1.03";
	const char *VersionDate = "November 20, 2011";

/* *** GLOBAL VARIABLES *** */

/* Mathematical Constants */
	const double Pi = 4.0 * atan(1.0);
	complex<double> ii = complex<double>(0.0,1.0);

/* *** Physical Parameters *** */

/* Electron Energy Parameter epsilon */
	double eps = 50.0;

/* Select Output Component:
		0 ... charge density, 
		1 ... current density (field direction),
		2 ... current density (radial direction),
		3 ... current density (azimuthal direction),
		4 ... current density (absolute value),
		5 ... current density (direction angle and magnitude in plane),
		6 ... charge velocity (current/density) 
*/ 
    int OutputType = 1;

/* *** Math Control Parameters *** */

/* Accuracy Goal */
	double Accu = 1.0e-10;

/* Maximum Landau Level */
	long MaxLevel = 10000;

/* *** Output Parameters (in scaled units) *** */

	double xMin = .0005;
	double xMax = 1.1;

	double zMin = -1.1;
	double zMax = 3.3;

/* Sampling Number */
  	int xSamples = 100;
	int zSamples = 400;

/* The output file name */
  	string BitmapFile = "At 50(quantum current test).bmp";

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
	double MaxDensity = 4.0 * eps;

/* Color Saturation Control */
	double Gamma = 0.5;
	double SaturationOffset = 1.0;

/* Integrated Current? */
	bool IntegratedCurrent = true;

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
		if (fabs(Value) > 1)
		{
			cerr << "ERROR:  Color Assignment Failure. \n\n";
			exit(1);
		}

	/* Blue Color Scheme for Positive Values, Red Color Scheme for Negative Values */
		if (Value >= 0)
		{
			*Red   = (unsigned char)floor(255.99 * (1.0 - sqrt(Value)));
			*Green = (unsigned char)floor(255.99 * (1.0 - Value));
			*Blue  = (unsigned char)floor(255.99 * (1.0 - Value * Value));
		}
		else
		{
			*Red   = (unsigned char)floor(255.99 * (1.0 - Value * Value));
			*Green = (unsigned char)floor(255.99 * (1.0 + Value));
			*Blue  = (unsigned char)floor(255.99 * (1.0 - sqrt(-Value)));
		}
	}

/* Determine Pixel Color from Complex Numerical Value [|z| <= 1] */
void AssignComplexColor(unsigned char *Red, unsigned char *Green, unsigned char *Blue, complex<double> Value)
{
	if (norm(Value) > 1)
	{
		cerr << "ERROR:  Color Assignment Failure. \n\n";
		exit(1);
	}

	/* Separate into Magnitude and Phase */
	double Magnitude = sqrt(norm(Value));
	double Phase = 0.0;
	if (Magnitude > 0.0)
	{
		Phase = arg(Value);
	}

	/* Translate Magnitude Into Brightness */
	double Brightness = 255.99 * Magnitude;

	/* Translate Phase Into Tint */
	double Aux = fabs(cos(0.5 * Phase));
	*Red   = (unsigned char)floor(Brightness * Aux);
	Aux = fabs(cos(0.5 * Phase - Pi / 3.0));
	*Green = (unsigned char)floor(Brightness * Aux);
	Aux = fabs(cos(0.5 * Phase + Pi / 3.0));
	*Blue  = (unsigned char)floor(Brightness * Aux);
}


/* *** MAIN PROGRAM ************************************************************** */

int _tmain(int argc, _TCHAR* argv[])
{

/* Opening Message */
	cerr << "MAGNGREEN.c++ - Energy Green Function in a Homogeneous Magnetic Field\n"
		 << "(Version " << VersionSpec << ", written " << VersionDate << ")\n\n";

/* *** STEP 1: CALCULATE A DENSITY PROFILE ***/

/* Variables */
	long i, j, n;
	double x, z, Aux, Aux2;
	complex<double> Sum, SumP, SumR, AuxC;

/* Start Message */
	cerr << "Calculating Density Profile: ";

/* Reserve Field for Density and Current Information */
	DensCurr *QuantumField = new DensCurr[xSamples * zSamples];

/* Reserve Array For 1D Green Functions */
	complex<double> *Green1DArray = new complex<double>[MaxLevel];
	complex<double> *Green1DPrimeArray = new complex<double>[MaxLevel];

/* Reserve and Calculate Array for Exponential Prefactors */
	double *Prefactor = new double[xSamples];

	for (j = 0; j < xSamples; ++j)
	{
		x = xMin + j * (xMax - xMin) / (xSamples - 1);

	/* Gaussian Envelope */
		Prefactor[j] = exp(-eps * x * x);
	}

/* Loop Through z Values */
	for (i = 0; i < zSamples; ++i)
	{
	/* Linear z Values */
		z = zMin + i * (zMax - zMin) / (zSamples - 1);

	/* Landau Level Counter */
		long Level = 0;

	/* Find Number of Open Channels */
		long Channel = (long)ceil((eps - 1.0) / 2.0);

	/* Green Function - Open Channels */
		while (Level < Channel)
		{
		/* Longitudinal Energy Parameter */
			Aux = sqrt(eps * (eps - (2 * Level + 1)));

		/* One-Dimensional Green Function in Field Direction, and its Derivative */
			Green1DPrimeArray[Level] = polar(1.0, 2.0 * Aux * fabs(z));
			Green1DArray[Level] = -ii * Green1DPrimeArray[Level] / Aux;

		/* Apply Sign Correction, If Necessary */
			if (z < 0) Green1DPrimeArray[Level] *= (-1);

		/* Increment Landau Level */
			++Level;
		}

	/* Green Function - Closed Channels */
   		do
		{
		/* Longitudinal Energy Parameter */
			Aux = sqrt(eps * ((2 * Level + 1) - eps));

		/* One-Dimensional Green Function in Field Direction, Derivative */
			Aux2 = exp(-2.0 * Aux * fabs(z));
			Green1DArray[Level] = - Aux2 / Aux;
			Green1DPrimeArray[Level] = Aux2;

		/* Apply Sign Correction, If Necessary */
			if (z < 0) Green1DPrimeArray[Level] *= (-1);

		/* Increment Landau Level */
			++Level;
		}
		while ((Aux2 > Accu) && (Level < MaxLevel));

	/* Indicate Progress and Convergence Failure */
		if (Level == MaxLevel)
			cerr << "!";
		else
			cerr << ".";

	/* Calculate Green Function, Derivative As Sums */
		for (j = 0; j < xSamples; ++j)
		{
			x = xMin + j * (xMax - xMin) / (xSamples - 1);

		/* Argument of Laguerre Polynomial */
			double Radial = 2.0 * eps * x * x;

		/* Initial Recursion Values For Laguerre Polynomials */
			double l0  = 1.0;
			double l1  = 1.0 - Radial;
			double l2;

		/* Same for Laguerre Derivative Function M_n = L_n + 2L_(n-1)^(1) */
			double m0  = l0;

		/* Calculate Green Function via Series */
			Sum = 0.0;
			SumP = 0.0;
			SumR = 0.0;

			for (n = 0; n < Level; ++n)
			{
			/* Add Term To Sum */
				Sum += l0 * Green1DArray[n];
				SumP += l0 * Green1DPrimeArray[n];
				SumR -= x * m0 * Green1DArray[n];

			/* Recursion Step For Laguerre Polynomials, Derivative Function */
				m0  += (l0 + l1);
				
				l2  = ((2 * n + 3 - Radial) * l1 - (n + 1) * l0) / (n + 2);
				l0  = l1;
				l1  = l2;
			}

		/* Apply Prefactors */
			Sum *= Prefactor[j];
			SumP *= Prefactor[j];
			SumR *= Prefactor[j];

		/* Integrated Density/Current, If Desired */
			Aux = 1.0;
			if (IntegratedCurrent == true)
			{
				Aux = 2.0 * Pi * x;
			}

		/* Pointer into Array */
			long Ptr = i * xSamples + j;
			
		/* Store Density, Current Density Vector */
			QuantumField[Ptr].Density          =  4.0 * Aux * eps * eps * norm(Sum);
			
			QuantumField[Ptr].CurrentField     =  4.0 * Aux * eps * imag(conj(Sum) * SumP);
			QuantumField[Ptr].CurrentRadial    =  4.0 * Aux * eps * eps * imag(conj(Sum) * SumR);
			QuantumField[Ptr].CurrentAzimuthal = -4.0 * Aux * eps * eps * x * norm(Sum);
		}
	}

/* Notify User */
	cerr << " done.\n\n";

/* Cleanup Arrays */
	delete[] Green1DArray;
	delete[] Green1DPrimeArray;
	delete[] Prefactor;

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
		/* Pointer into Array */
			long Ptr = i * xSamples + j;

		/* Select Desired Output */
			bool PlotReal = true;
			
			switch (OutputType)
			{
			/* Charge Density */
			case 0:
				Aux = QuantumField[Ptr].Density;
				break;
			/* Current Density (field direction) */
			case 1:
				Aux = QuantumField[Ptr].CurrentField;
				break;
			/* Current Density (radial direction) */
			case 2:
				Aux = QuantumField[Ptr].CurrentRadial;
				break;
			/* Current Density (azimuthal direction) */
			case 3:
				Aux = QuantumField[Ptr].CurrentAzimuthal;
				break;
			/* Current Density (absolute value) */
			case 4:
				Aux = sqrt(QuantumField[Ptr].CurrentField * QuantumField[Ptr].CurrentField 
					+ QuantumField[Ptr].CurrentRadial * QuantumField[Ptr].CurrentRadial 
					+ QuantumField[Ptr].CurrentAzimuthal * QuantumField[Ptr].CurrentAzimuthal);
				break;
			/* Current Density (direction angle map) */
			case 5:
				PlotReal = false;
				AuxC = QuantumField[Ptr].CurrentRadial + ii * QuantumField[Ptr].CurrentField;
				break;
			/* "Charge Velocity" (absolute value of current / density) */
			case 6:
				Aux = sqrt(QuantumField[Ptr].CurrentField * QuantumField[Ptr].CurrentField 
					+ QuantumField[Ptr].CurrentRadial * QuantumField[Ptr].CurrentRadial 
					+ QuantumField[Ptr].CurrentAzimuthal * QuantumField[Ptr].CurrentAzimuthal) / QuantumField[Ptr].Density;
				break;
			default:
				cerr << "ERROR:  Invalid Output Type " << OutputType << "\n\n";
				exit(0);
			}

		/* Select Plot Type */
			if (PlotReal == true)
			{
			/* Real Plots */

			/* Calculate Relative Density */
				Aux  /= MaxDensity;	

			/* Exposure Limits */
				if (Aux > 1.0) Aux = 1.0;
				if (Aux < -1.0) Aux = -1.0;

			/* Sort Out Unreliable Results */
				if (_isnan(Aux) != 0.0)
				{
					cerr << "Warning:  Calculation failed.\n";

					Red = (unsigned char)0;
					Green = (unsigned char)255;
					Blue = (unsigned char)0;
				}
			else
				{
				/* Adjust Color Intensity */
					Saturation = pow(fabs(Aux), Gamma);
					if (Aux < 0) Saturation *= (-1);

				/* Assign Coloring */
					AssignColor(&Red, &Green, &Blue, Saturation);
				}
			}
			else
			{
			/* Complex Plot */

			/* Calculate Relative Density */
				AuxC /= MaxDensity;

			/* Apply Saturation Offset for Complex Graphs */
				AuxC *= SaturationOffset;

			/* Adjust Color Intensity */
				Aux = sqrt(norm(AuxC));
				Saturation = pow(Aux, Gamma);
				AuxC = AuxC * Saturation / Aux;
				
			/* Exposure Limits (Accu = small parameter to prevent numerical overflow) */
				while (norm(AuxC) > 1.0)
				{
					AuxC = AuxC / sqrt(norm(AuxC) + Accu);
				}

			/* Sort Out Unreliable Results */
				if ((_isnan(real(AuxC)) != 0.0) || (_isnan(imag(AuxC)) != 0.0))
				{
					cerr << "Warning:  Calculation failed.\n";

					Red = (unsigned char)255;
					Green = (unsigned char)255;
					Blue = (unsigned char)255;
				}
				else
				{
				/* Assign Coloring */
					AssignComplexColor(&Red, &Green, &Blue, AuxC);
				}
			}

		/* Paint Pixel (Note: Progress top-to-bottom) */
			PSet(ImageTemplate, BMPWidth, BMPHeight, j, zSamples - i - 1, Red, Green, Blue);
		}

	/* Progress Marker */
		cerr << ".";
	}

/* Clean Up Fields */
	delete[] QuantumField;

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