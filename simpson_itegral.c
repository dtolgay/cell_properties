#include <stdio.h>


double simpson_integral(double *y, double *x, const int size)
{

    /*This function is created to take the simpson integral of the emissivity values in the radiationField_J.dat file.
    If the number of intervals is even then the Simpsons 1/3 rule is used and if the number of intervals is odd Simpson's 1/3 rule and Simpson's 3/8 rule is combined to calculate the integral.

    Parameters:
    ------------------- 
    y: array 

    x: array

    size: int 

    Return:
    ------------------- 
    result: double 
    This is the integrated value.
    

    */

    // Initialize the result 
    double result = 0;

    // If the size of the array is even, then the number of intervals that the integration will be calculated is an odd number.
    // If number of interval is an odd number one cannot use the Simpsons 1/3 rule directly. 
    // The code below first calculates the integral for the size - 3 term. Then the remaining 4 terms y[size-4], y[size-3], y[size-2], y[size-1] are being calculated using the Simpson's 3/8 rule.
    if (size % 2 == 0)
    {
        double const step_size = x[1] - x[0];

        // This is the first part of the Simpson's 1/3 integral
        int i = 0;
        for ( ; i < size - 3; i++)
        {

            if (i == 0 || i == size - 4)  // Only for the start of this portion
            {
                result += y[i];
            }
            else if (i % 2 == 1)  // Odd index in this portion
            {
                result += 4 * y[i];
            }
            else  // Even index in this portion
            {
                result += 2 * y[i];
            }
        }

        // Multiply the result of Simpson's 1/3 part by its coefficient
        result = result * step_size / 3.0;

        // This is the remaining integration for Simpson's 3/8 rule
        double simpson_3_8 = (y[i-1] + 3 * y[i] + 3 * y[i+1] + y[i+2]) * step_size * 3.0 / 8.0; // Starting with y[i-1] for the overlapping point
        result += simpson_3_8;
        return result; 
    }

    
    // If the size of the array is odd, then it means that there are even number of intervals there and Simpson's 1/3 rule can be used throughly to calculate the integral.
    if (size % 2 == 1)
    {
        // Integrate according to the Simpson's 1/3 rule

        double const step_size = x[1] - x[0];  // Assuming uniform step size

        for (int i = 0; i < size; i++)
        {
            if (i == 0 || i == size-1)
            {
                result += y[i];
            }
            else if (i % 2 == 1)
            {
                result += 4 * y[i];
            }
            else
            {
                result += 2 * y[i];
            }
        }  
        result *= step_size / 3.0;
        return result; 
    }
}

int main()
{
	double y[99];	
	double x[99];

	for (int i = 0; i < 99; i++)
	{
		x[i] = 100 + i;
		y[i] = 1/x[i];
	}

	printf("x[0] = %f ---- x[98] = %f\n", x[0], x[98]);


	double result = simpson_integral(y, x, 99);

	printf("%f\n", result);

	return 0;

}