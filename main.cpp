#include <cmath>
#include <iomanip>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "polyfit.c"
#include "testADC.h"
#include "stdlib.h"
#include "memory.h"

using namespace std;

const int ADCREVS          = 5;
const unsigned int         RANGE = 200;
const unsigned int         TRIM  = 200;
const unsigned int         SKIP_RANGE  = 50;

enum
{
    LINEAR
};

CADC* pADC;
// cross correlation in 2D vector by inputA(row), inputB(other row)
// range = 200
// trim = 200

int FindIndex(std::vector<double> alist, double item)
{
    if ( (item < *(alist.begin())) || (item > *(alist.end() - 1)) )
    {
        return -1;
    }
    else if (item == *(alist.begin()))
    {
        return 0;
    }
    else if (item == *(alist.end() - 1))
    {
        return (alist.size() - 1);
    }
    else
    {
        int Min_idx = 0;
        int Max_idx = (alist.size() - 1);
        int Mid_idx = 0;
        while ((Max_idx - Min_idx) > 1)
        {
            Mid_idx = int((Max_idx + Min_idx) / 2);
            if (item < *(alist.begin() + Mid_idx))
            {
                Max_idx = Mid_idx;
            }
            else
            {
                Min_idx = Mid_idx;
            }
        }
        return Min_idx;
    }
}

int vXcorr(const std::vector< std::vector<double> > &pVec, unsigned int inputA, unsigned int inputB, std::vector<double> &pOutputA, std::vector<double> &pOutputB)
{
   int i, ri = 0;
   unsigned int j;
   double xySigma, xSigma, ySigma, xPowSigma, yPowSigma, maxCorr;
   int idx;
   std::vector <double> correlation;

   maxCorr = 0;

   for (i = -30; i < static_cast<int>(RANGE - 30); i++)
   {
      idx = 0;
      xySigma = 0;
      xSigma = 0;
      ySigma = 0;
      xPowSigma = 0;
      yPowSigma = 0;
      for (j = 0 + SKIP_RANGE; j < RANGE+SKIP_RANGE; j++)
      {
         xySigma += (pVec[inputA][j] *  pVec[inputB][i+j]);
         xSigma += pVec[inputA][j];
         ySigma += pVec[inputB][i+j];
         xPowSigma += pow(pVec[inputA][j], 2);
         yPowSigma += pow(pVec[inputB][i+j], 2);

         idx++;
      }

      correlation.push_back(((idx * xySigma) - (xSigma * ySigma)) / sqrt(((idx * xPowSigma) - pow(xSigma, 2)) * ((idx * yPowSigma) - pow(ySigma, 2))));

      if (maxCorr < correlation.back())
      {
         maxCorr = correlation.back();
         ri = i;
      }
   }// for

   for (unsigned int k = TRIM, idx = 0; k < TRIM + pVec[inputA].size(); k++, idx++)
   {
      j = k + ri;

      if (inputB == 1)
      {
         pOutputA.push_back(pVec[inputA][k]);
      }
      pOutputB.push_back(pVec[inputB][j]);
   }

   return ri;
}


void cppInsert()
{
   pADC->init_ADC_vector(ADCREVS);
}

void cppDelete()
{
   delete pADC;
}

void testVector( std::vector<double> &pTest )
{
      std::cout << pTest[1] << " ";

      pTest[1] = 18  ;
      pTest.push_back(11);
}

std::vector<double> averageZeroToBeforelast(std::vector< std::vector<double> > &pVector, int option = 1)
{
    std::vector<double> pReturn;
    double sum;
    unsigned int rows;
    unsigned int cols = pVector.at(0).size();

    if (option == 1)
    {
       rows = pVector.size()-1; // except last one.
    }
    else
    {
       rows = pVector.size();
    }

    for (unsigned int i = 0; i<cols; i++)
    {
      sum = 0;
       for (unsigned int j = 0; j<rows; j++)
       {
          sum += pVector[j][i];
       }

       pReturn.push_back(sum/static_cast<double>(rows));
    }

    return pReturn;
}

#define INVALID_DATA     1818181818.0

std::vector<double> interp1(std::vector<double> X, std::vector<double> V, std::vector<double> Xq, int rMode)
{
    std::vector<double> Rst;
    if (rMode == LINEAR)
    {
        int i = 0;
        for (std::vector<double>::iterator it_Xq = Xq.begin(); it_Xq != Xq.end(); ++it_Xq)
        {
            if ( (*it_Xq < *(X.begin())) || (*it_Xq > *(X.end() - 1)) )
            {
                Rst.push_back(INVALID_DATA);
                i++;
            }
            else
            {
                int j = FindIndex(X, Xq[i]);

                double y0 = V[j];
                double y1 = V[j + 1];
                double x = Xq[i++];
                double x0 = X[j];
                double x1 = X[j + 1];

                double result = y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
                Rst.push_back(result);
            }
        }
    }
    return Rst;
}

//////////////////////////////////////////////////////////////////////
//////////////////// post processing /////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


   const unsigned int         LENGTH_PREAMPLE   = 500;
   const unsigned int         LENGTH_PROCESS    = 30000;
   const unsigned int         LENGTH_PRBS       = 255;
   const unsigned int         OVERSAMPLING      = 2;
   const unsigned int         LENGTH_PRBS_PERIOD= round(LENGTH_PRBS*OVERSAMPLING*1.1);


int main()
{

   string line;
   ifstream myfile( "After_Brickwall.txt" );

   std::vector< std::vector <double> > vADC;
   std::vector< std::vector <double> > xCorr;
   std::vector< std::vector <double> > xAligned;

   vADC.resize( READ_ITR, std::vector<double>() );
   xCorr.resize(READ_ITR, std::vector<double>() );
   xAligned.resize(READ_ITR, std::vector<double>() );

   std::vector<double> vAvg;
   std::vector<double> vTemp;
   std::vector<double> vX;
   std::vector<double> vXQ;
   std::vector<double> vShift;

   // for input from file
   unsigned int   nRow     = 0;
   unsigned int   curCol   = 0;
   double         temp;
   double   *CoefsReverse= NULL;    // polyfit


   CoefsReverse = (double *)malloc((3)*sizeof (double));

   // input xcorred adc
   if (myfile)  // same as: if (myfile.good())
   {
      while (getline( myfile, line ))  // same as: while (getline( myfile, line ).good())
      {
         // line is a row (1 row = 50cols)
         // matlab raw data file = 1 col is 1 sector read adc   first row = first adc data of 50 sector read adc
         if (nRow < POSTPRCESSSIZE)
         {
            //cout << line << endl;
            nRow++;

            istringstream row(line);

            curCol = 0;
            while (row)
            {
               string cols;

               if (!getline(row, cols, ',')) break;

               if (curCol < READ_ITR)
               {
                  temp = strtod(cols.c_str(), NULL);
                  //xCorr[curCol].push_back(temp);
                  vADC[curCol].push_back(temp);

                  curCol++;
               }
               else
               {
                  break;
               }
            }
         }
         else
         {
            break;
         }

      }

      cout << "rows= " << nRow << endl;
      cout << "cols= " << curCol << endl;
      myfile.close();
   }
   else cout << "error -_-\n";

   // Xcorr
   unsigned int adcRows = vADC.size();
   unsigned int adcCols = vADC[0].size();

   cout << "rows= " << adcRows << endl;
   cout << "cols= " << adcCols << endl;

   for (unsigned int i=1; i<adcRows; i++)
   {
      if (i==1)
      {
         vXcorr(vADC, 0, 1, xCorr[0], xCorr[1]);
      }
      else
      {
         vXcorr(vADC, 0, i, xCorr[READ_ITR], xCorr[i]);
      }
   }

   // save output file
   unsigned int rowLength = xCorr.size();
   unsigned int colLength = xCorr[0].size();

   cout << "rowLength= " << rowLength<< endl;
   cout << "colLength= " << colLength<< endl;


   ofstream rstfile;
   rstfile.open ("FineAlignedSignal.txt");
   std::cout << "Create output file FineAlignedSignal.txt" << std::endl;


   for (unsigned int i=0; i<colLength; i++)
   {
      for (unsigned int j=0; j<rowLength; j++)
      {
         if (j == rowLength-1)
         {
            rstfile << setprecision(12) << xCorr[j][i];
         }
         else
         {
            rstfile << setprecision(12) << xCorr[j][i] << ", ";
         }

      }
      rstfile << std::endl;
   }


   rstfile.close();

   // average all xcorr adc
    vAvg = averageZeroToBeforelast(xCorr, FULLAVERAGE); // real code
   for (unsigned int i=0; i<POSTPRCESSSIZE; i++)
   {
      //vAvg.push_back(yt[i]);    // only test vAvg
      vX.push_back(i);
   }

   // loop read itr
   // first loop for generate signal alignment coefficient
   vShift.clear();

   for (unsigned int i=0; i<READ_ITR; i++)  // read loop
   {
      std::vector<double> vInterpl1;

      memset(Ypoints, 0, LENGTH_COEF);

      for (unsigned int j=0; j<LENGTH_COEF; j++)  // shift coefficient
      {
         vXQ.clear();
         for (unsigned int k=0; k<POSTPRCESSSIZE; k++)
         {
            vXQ.push_back(k + shiftCoef[j]);
         }

         vInterpl1.clear();
         vInterpl1 = interp1(vX, xCorr[i], vXQ, LINEAR);

         double delta = 0;
         vTemp.clear();
         for (unsigned int idxAvg=0; idxAvg<vAvg.size(); idxAvg++)
         {
            if (vInterpl1[idxAvg] == INVALID_DATA)
            {
               delta = 0;
            }
            else
            {
               delta = vAvg[idxAvg] - vInterpl1[idxAvg];
            }

            vTemp.push_back(delta);
         }

         // .2
         for (unsigned int idxTemp=0; idxTemp<vTemp.size(); idxTemp++)
         {
            vTemp[idxTemp] = vTemp[idxTemp] * vTemp[idxTemp];
         }

         // sum
         double sumFx = 0;
         for (unsigned int idxSum=1; idxSum<vTemp.size()-1; idxSum++)
         {
            sumFx += vTemp[idxSum];
         }

         // end of fx-loop
         Ypoints[j] = sumFx;
         std::cout << std::setprecision(12) << Ypoints[j] << std::endl ;
      }

      // polyfit !
      poly_fit(2, CoefsReverse, shiftCoef, Ypoints, LENGTH_COEF, NULL, NULL);

      double tempShift = -CoefsReverse[1] / (2.0 * CoefsReverse[2]);

      vShift.push_back(tempShift);
      cout << "Coeffs =" << std::setprecision(12) << CoefsReverse[0] << " " << CoefsReverse[1] << " " << CoefsReverse[2] << " " << " Shift = " << tempShift << std::endl;
   }

   // find signal alignment
   std::vector<double> vInterpl2;

   for (unsigned int i=0; i<READ_ITR; i++)
   {
      std::cout << "shift = " <<  i << " " <<std::setprecision(12) << vShift[i] << std::endl ;
   }

   for (unsigned int i=0; i<READ_ITR; i++)  // read loop
   {
      vInterpl2.clear();
      vXQ.clear();
      for (unsigned int k=0; k<POSTPRCESSSIZE; k++)
      {
         vXQ.push_back(k + vShift[i]);
      }

      vInterpl2 = interp1(vX, xCorr[i], vXQ, LINEAR);
      vInterpl2.erase(vInterpl2.begin());
      vInterpl2.erase(vInterpl2.end()-1);
      xAligned[i] = vInterpl2;

      std::cout << "Fine alignment " << i << " interpolation done" << std::endl;
   }

   // save output file
   for (unsigned int i=0; i<colLength; i++)
   {
      for (unsigned int j=0; j<rowLength; j++)
      {
         if (j == rowLength-1)
         {
            rstfile << setprecision(12) << xAligned[j][i];
         }
         else
         {
            rstfile << setprecision(12) << xAligned[j][i] << ", ";
         }

      }
      rstfile << std::endl;
   }

   rstfile.close();

   std::cout << "Processing End" << std::endl;
   free(CoefsReverse);

   return 0;
}
