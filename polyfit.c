#ifndef POLYFIT_CPP_INCLUDED
#define POLYFIT_CPP_INCLUDED


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>


#define POLY_FIT_ERROR_LIMIT_EXCEEDED     48410    /* Polynomial Fit Error Limit Exceeded */
#define OUT_OF_MEMORY                     10532   /* Tester Comm Mgr-Can't Allocate Memory */                  // SF2/3 various
#define DBL_EPSILON         2.2204460492503131e-16
#define FLT_EPSILON         1.19209290e-7F
#define LDBL_EPSILON        2.2204460492503131e-16L

#define f0pz     0.0     /* double */
#define f1pz     1.0     /* double */
#define f6E8     6.0E8   /* double */
#define f11Em2   0.110   /* double */

#define D0_0           0.0   // double
#define D0_5           0.5   // double
#define D1_0           1.0   // double
#define D2_0           2.0   // double



typedef short              sint16;
typedef unsigned short uint16;
typedef unsigned int uint32;

void free_double_matrix(double **matrix);
void LU_back_substitute(double **matrix, uint16 nn, sint16 *perm, double solution[]);
sint16 LU_decompose(double **matrix, uint16 nn, sint16 *perm, double *parity);
void form_dependent_knowns_double_vector(double *pvector, sint16 vector_dim, double *independent_data, double *dependent_data, sint16 nbr_data_points);
sint16 poly_fit(sint16 poly_order, double *solution, double *independent_data, double *dependent_data,sint16 nbr_data_points, double *Rsquared, double *std_dev);
void LU_back_sub_dbl(double **matrix, uint16 nn, sint16 *perm, double solution[]);


void free_double_matrix(double **matrix)
{
   if (matrix)
   {
      if (matrix[0])
      {
         free(matrix[0]);  // First free matrix space.
      }
      free(matrix);     // Then free row pointer space.
   }
}


void LU_back_substitute(double **matrix, uint16 nn, sint16 *perm, double solution[])
{
   sint16 row;      // matrix row
   sint16 col;      // matrix column
   sint16 ii = -1;  // interchange row
   sint16 ip;       // interchange permuation
   double sum;       // sum of the left hand sides
   double tempf;

   // Perform forward substitution
   for (row = 0; row < nn; row++)
   {
      ip = perm[row];
      sum = solution[ip];
      solution[ip] = solution[row];
      if (ii >= 0)
      {
         for (col = ii; col <= row - 1; col++)
         {
            tempf = matrix[row][col] * solution[col];
            sum = sum - tempf;
         }
      }
      else if (sum)
      {
         ii = row;
      }
      solution[row] = sum;
   }

   // Perform back substitution
   for (row = nn - 1; row >= 0; row--)
   {
      sum = solution[row];
      for (col = row + 1; col < nn; col++)
      {
         tempf = matrix[row][col] * solution[col];
         sum = sum - tempf;
      }
      solution[row] = sum / matrix[row][row];
   }
}




sint16 LU_decompose(double **matrix, uint16 nn, sint16 *perm, double *parity)
{
   sint16 row;          // matrix row
   sint16 col;          // matrix column
   sint16 kk;
   sint16 max_row = 0;  // Initialize to prevent compiler warning.
   double mabs;          // Maximum absolute value of current row;
   double dum;
   double sum;
   double temp;
   double *scale;        // stores the implicit scaling of each row.
   double tempf;
   sint16 rtrn_stat = 0;

   scale = (double *)malloc(nn * sizeof(*scale));
   if (!scale)
   {
      rtrn_stat = OUT_OF_MEMORY;
      goto LU_decompose_exit;
   }

   *parity = 1.0f;                  // Set parity to even to indicate no row interchanges have occurred.

   for (row = 0; row < nn; row++)   // Loop over rows to get the implicit scaling information.
   {
      mabs = 0.0f;
      for (col = 0; col < nn; col++)
      {
         if ((temp = fabs(matrix[row][col])) > mabs)
            mabs = (double)temp;
      }
      if (mabs == 0.0f)
      {
         rtrn_stat = (sint16)POLY_FIT_ERROR_LIMIT_EXCEEDED;   // Was SINGULAR_MATRIX, replace with valid codes.h return code
         goto LU_decompose_exit;
      }
      scale[row] = 1.0f / mabs;         // Implied scale factors by row; mabs = Maximum absolute value in the row.
   }

   for (col = 0; col < nn; col++)
   {
      for (row = 0; row < col; row++)
      {
         sum = matrix[row][col];
         for (kk = 0; kk < row; kk++)
         {
            tempf = matrix[row][kk] * matrix[kk][col];
            sum -= tempf;
         }
         matrix[row][col] = sum;
      }

      mabs = 0.0f;
      for (row = col; row < nn; row++)
      {
         sum = matrix[row][col];
         for (kk = 0; kk < col; kk++)
         {
            tempf = matrix[row][kk] * matrix[kk][col];
            sum -= tempf;
         }
         matrix[row][col] = sum;
         if ((dum = scale[row] * (double)fabs(sum)) >= mabs)
         {
            mabs = dum;
            max_row = row;
         }
      }

      if (col != max_row)
      {
         for (kk = 0; kk < nn; kk++)
         {
            dum = matrix[max_row][kk];
            matrix[max_row][kk] = matrix[col][kk];
            matrix[col][kk] = dum;
         }
         *parity = -(*parity);
         scale[max_row] = scale[col];
      }

      perm[col] = max_row;

      if (matrix[col][col] == 0.0f)
         matrix[col][col] = FLT_EPSILON;     // Substitute the smallest effective value to prevent divide by zero.

      if (col != (nn - 1))
      {
         dum = 1.0f / (matrix[col][col]);
         for (row = col + 1; row < nn; row++)
            matrix[row][col] *= dum;
      }
   }
 LU_decompose_exit:
   if (scale)
      free(scale);

   return(rtrn_stat);
}


void form_dependent_knowns_double_vector(double *pvector, sint16 vector_dim, double *independent_data,
                                         double *dependent_data, sint16 nbr_data_points)
{
   sint16 row;       // vector row
   sint16 ii;
   double tempf;

   // Form column vector of dependent knowns.
   for (row = 0; row < vector_dim; row++)
   {
      pvector[row] = f0pz;
      for (ii = 0; ii < nbr_data_points; ii++)
      {
         tempf = (double)pow(independent_data[ii], row) * dependent_data[ii];
         pvector[row] += tempf;
      }
   }
}




void form_independent_knowns_double_matrix(double **matrix, sint16 matrix_dim, double *independent_data,
                                          sint16 nbr_data_points)
{
   sint16 row;       // matrix row
   sint16 col;       // matrix column
   sint16 ii;
   double tempf;

   // Form matrix of independent known values.
   for (row = 0; row < matrix_dim; row++)
   {
      for (col = 0; col < matrix_dim; col++)
      {
         matrix[row][col] = 0.0f;
         for (ii = 0; ii < nbr_data_points; ii++)
         {
            tempf = (double)pow(independent_data[ii], (row + col));
            matrix[row][col] += tempf;
         }
      }
   }
}



double  d_abs( double input)
{
    double output;
    if (input < 0)
        output = -input;
    else
        output = input;

    return output;
}

double **alloc_double_matrix(uint16 rows, uint16 columns)
{
   double **matrix = NULL;
   sint16 row;

   // Allocate space for row pointers.
   matrix = (double **)malloc(rows * sizeof(double *));
   if (matrix == NULL)
   {
      return NULL;
   }

   // Allocate matrix space. Assign resulting pointer to 1st row.
   matrix[0] = (double *)malloc((uint32)rows * (uint32)columns * sizeof(double));
   if ( matrix[0] == NULL )
   {
      free(matrix);
      return NULL;
   }

   // Set matrix row pointers.  Start at 1 since we assign 0 during the malloc
   for ( row = 1; row < rows; row++ )
   {
      matrix[row] = matrix[0] + (row * columns);
   }
   return (matrix);
}



sint16 LU_decompose_dbl(double **matrix, uint16 nn, sint16 *perm, double *parity)
{
   sint16 row;          // matrix row
   sint16 col;          // matrix column
   sint16 kk;
   sint16 max_row = 0;  // Initialize to prevent compiler warning.
   double mabs;          // Maximum absolute value of current row;
   double dum;
   double sum;
   double temp;
   double *scale;        // stores the implicit scaling of each row.
   double tempf;
   sint16 rtrn_stat = 0;

   scale = (double *)malloc(nn * sizeof(*scale));
   if (!scale)
   {

      goto LU_decompose_exit;
   }

   *parity = 1.0;                  // Set parity to even to indicate no row interchanges have occurred.

   for (row = 0; row < nn; row++)   // Loop over rows to get the implicit scaling information.
   {
      mabs = 0.0;
      for (col = 0; col < nn; col++)
      {
         if ((temp = d_abs(matrix[row][col])) > mabs)
            mabs = (double)temp;
      }
      if (mabs == 0.0)
      {
         rtrn_stat = (sint16)POLY_FIT_ERROR_LIMIT_EXCEEDED;   // Was SINGULAR_MATRIX, replace with valid codes.h return code
         goto LU_decompose_exit;
      }
      scale[row] = 1.0 / mabs;         // Implied scale factors by row; mabs = Maximum absolute value in the row.
   }

   for (col = 0; col < nn; col++)
   {
      for (row = 0; row < col; row++)
      {
         sum = matrix[row][col];
         for (kk = 0; kk < row; kk++)
         {
            tempf = matrix[row][kk] * matrix[kk][col];
            sum -= tempf;
         }
         matrix[row][col] = sum;
      }

      mabs = 0.0;
      for (row = col; row < nn; row++)
      {
         sum = matrix[row][col];
         for (kk = 0; kk < col; kk++)
         {
            tempf = matrix[row][kk] * matrix[kk][col];
            sum -= tempf;
         }
         matrix[row][col] = sum;
         if ((dum = scale[row] * d_abs(sum)) >= mabs)
         {
            mabs = dum;
            max_row = row;
         }
      }

      if (col != max_row)
      {
         for (kk = 0; kk < nn; kk++)
         {
            dum = matrix[max_row][kk];
            matrix[max_row][kk] = matrix[col][kk];
            matrix[col][kk] = dum;
         }
         *parity = -(*parity);
         scale[max_row] = scale[col];
      }

      perm[col] = max_row;

      if (matrix[col][col] == 0.0)
         matrix[col][col] = DBL_EPSILON;     // Substitute the smallest effective value to prevent divide by zero.

      if (col != (nn - 1))
      {
         dum = 1.0 / (matrix[col][col]);
         for (row = col + 1; row < nn; row++)
            matrix[row][col] *= dum;
      }
   }
 LU_decompose_exit:
   if (scale)
      free(scale);

   return(rtrn_stat);
}


void LU_back_sub_dbl(double **matrix, uint16 nn, sint16 *perm, double solution[])
{
   sint16 row;      // matrix row
   sint16 col;      // matrix column
   sint16 ii = -1;  // interchange row
   sint16 ip;       // interchange permuation
   double sum;       // sum of the left hand sides
   double tempf;

   // Perform forward substitution
   for (row = 0; row < nn; row++)
   {
      ip = perm[row];
      sum = solution[ip];
      solution[ip] = solution[row];
      if (ii >= 0)
      {
         for (col = ii; col <= row - 1; col++)
         {
            tempf = matrix[row][col] * solution[col];
            sum = sum - tempf;
         }
      }
      else if (sum)
      {
         ii = row;
      }
      solution[row] = sum;
   }

   // Perform back substitution
   for (row = nn - 1; row >= 0; row--)
   {
      sum = solution[row];
      for (col = row + 1; col < nn; col++)
      {
         tempf = matrix[row][col] * solution[col];
         sum = sum - tempf;
      }
      solution[row] = sum / matrix[row][row];
   }
}


sint16 poly_fit(sint16 poly_order, double *solution, double *independent_data, double *dependent_data,
               sint16 nbr_data_points, double *Rsquared, double *std_dev)
{
   sint16 return_stat = 0;
   double **matrix = NULL;
   sint16 *perm = NULL;
   double parity;
   sint16 matrix_dim = poly_order + 1;
   sint16 ii;
   sint16 jj;
   double fitted_dependent_value;
   double error;
   double sum_squared_error;
   double mean;
   double mean_delta;
   double sum_squared_mean_delta;
   double tempf;

   matrix = alloc_double_matrix(matrix_dim, matrix_dim);
   if (!matrix)
   {
      return_stat = OUT_OF_MEMORY;
      goto poly_fit_exit;
   }
   perm = (sint16 *)malloc(matrix_dim * sizeof(*perm));
   if (!perm)
   {
      return_stat = OUT_OF_MEMORY;
      goto poly_fit_exit;
   }

   form_independent_knowns_double_matrix(matrix, matrix_dim, independent_data, nbr_data_points);
   form_dependent_knowns_double_vector(solution, matrix_dim, independent_data, dependent_data, nbr_data_points);

   // Decompose normalized cylinders matrix to L*U format.
   return_stat = LU_decompose(matrix, matrix_dim, perm, &parity);

   // Perform substitutions for the sample set.
   // NOTE: This leaves the polynomial coefficients in "reverse" order.
   if (!return_stat)
      LU_back_substitute(matrix, matrix_dim, perm, solution);

   // Initialize variables outside the conditional, to avoid uninitialize variable compiler warning.
   sum_squared_error = 0.0f;
   sum_squared_mean_delta = 0.0f;

   if (Rsquared || std_dev)    // Calculate the "Coefficient of Determination" OR the standard deviation.
   {
      mean = 0.0f;

      // Calculate measured dependent mean.
      for (ii = 0; ii < nbr_data_points; ii++)
      {
         mean += dependent_data[ii];
      }
      mean /= nbr_data_points;

      for (ii = 0; ii < nbr_data_points; ii++)
      {
         fitted_dependent_value = solution[0];
         for (jj = 1; jj < matrix_dim; jj++)
         {
            tempf = solution[jj] * (double)pow(independent_data[ii], jj);
            fitted_dependent_value += tempf;
         }
         mean_delta = dependent_data[ii] - mean;
         sum_squared_mean_delta += mean_delta * mean_delta;
         error = dependent_data[ii] - fitted_dependent_value;
         sum_squared_error += error * error;
      }
   }

   if (Rsquared)
   {
      if (sum_squared_mean_delta == 0.0)
      {
         *Rsquared = 1.0;
      }
      else
      {
         *Rsquared = 1.0f - (sum_squared_error / sum_squared_mean_delta);
      }
   }


   if (std_dev)
   {
      *std_dev = (double)sqrt(sum_squared_error / (nbr_data_points - (poly_order+1)));
   }

  poly_fit_exit:
   if (perm)
      free(perm);
   if (matrix)
      free_double_matrix(matrix);
   return(return_stat);
}



#endif // POLYFIT_CPP_INCLUDED


