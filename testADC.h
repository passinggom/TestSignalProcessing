#ifndef TESTADC_H_INCLUDED
#define TESTADC_H_INCLUDED

#include <iostream>
#include <vector>

const unsigned int         READ_ITR          = 50;
const unsigned int         LENGTH_COEF       = 9;
//const unsigned int         POSTPRCESSSIZE    = 27970;
const unsigned int         POSTPRCESSSIZE    = 29001;
const int                  FULLAVERAGE       = 0;

double  shiftCoef[LENGTH_COEF] = {-1, -0.75, -0.5, -0.25, 0,  0.25, 0.5, 0.75, 1};
double  Ypoints[LENGTH_COEF];

// the actual class
class CADC
{
private:
   std::vector< std::vector<int> > vTest;        // 2D vector
   unsigned int adcRevs;                        // matrix row size
   unsigned int initOrNot;                      // Before use CADC class, you should be call init_ADC_vector;

public:
    std::vector< std::vector<int> > vAdc;        // 2D vector
    void             init_ADC_vector   (unsigned int revs);
    void             insertADC         (unsigned int revs, int value);
    unsigned int     checkADCInit      (void);
    void             setADCInit        (unsigned int flag);

    unsigned int     getVectorRowSize  (void);     // all row is same size
    unsigned int     getVectorColSize  (void);     // all col is same size

    ~CADC();
    CADC();
};

CADC::CADC()
{

}

CADC::~CADC()
{
    std::vector< std::vector<int> >().swap(vAdc);    // release memory
    std::vector< std::vector<int> >().swap(vTest);    // release memory
}

unsigned int CADC::getVectorRowSize (void)
{
   return vAdc.size();
}

unsigned int CADC::getVectorColSize (void)
{
   return vAdc[0].size();
}


void CADC::init_ADC_vector (unsigned int revs)
{
   adcRevs = revs;
   setADCInit(1);
   vAdc.resize( adcRevs , std::vector<int>( ) );
}

void CADC::insertADC (unsigned int revs, int value)
{
   vAdc.at(revs).push_back(value);
}

unsigned int CADC::checkADCInit (void)
{
    return initOrNot;
}

void CADC::setADCInit (unsigned int flag)
{
    initOrNot = flag;
}


#endif // TESTADC_H_INCLUDED
