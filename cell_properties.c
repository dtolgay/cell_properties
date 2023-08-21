#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Define the number of columns in the radiationField_J.dat and radiationField_wavelenghts.dat 
#define NUM_COLUMNS_J 11
#define NUM_COLUMNS_WAVELENGTHS 4

// Define constants
const double c = 3e8;               // Speed of light:          m/s
const double u_HAB = 5.29e-15;      // Habing energy density:   J/m^3
const double G0 = 1.69;             // Habing parameter:        Unitless
const double pi = 3.14159265359;    // pi constant

int countDataRows(const char *filename) {
    int count = 0;
    char line[1024];
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '#') { // Exclude lines starting with '#'
            count++;
        }
    }
    
    fclose(fp);
    return count;
}

int readRadiationField_J_dat_File(const char *filename, double **radiationField)
{

    int rowNum = 0;

    // Open the file
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    char line[1024];  // buffer for each line

    while (fgets(line, sizeof(line), fp)) {
        // Skip lines starting with '#'
        if (line[0] == '#') {
            continue;
        }

        // Parse the data
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                   &radiationField[rowNum][0], &radiationField[rowNum][1], 
                   &radiationField[rowNum][2], &radiationField[rowNum][3], 
                   &radiationField[rowNum][4], &radiationField[rowNum][5], 
                   &radiationField[rowNum][6], &radiationField[rowNum][7], 
                   &radiationField[rowNum][8], &radiationField[rowNum][9], 
                   &radiationField[rowNum][10]) == 11) {
            rowNum++;
        }
    }

    fclose(fp);

    return 0;
}

int readWWavelenghts(const char *filename, double **wavelengths)
{
    int rowNum = 0; 

    //Open the file 
    FILE *fp = fopen(filename, "r");
    // Check there is a memory that file pointer points to
    if (fp == NULL)
    {
        perror("Unable to open file!"); 
        exit(1);
    }

    // Cheate a buffer for each line 
    char line[1024]; 

    // Read the file 
    while (fgets(line, sizeof(line), fp))
    {
        // Skip the header lines 
        if (line[0] == '#')
        {
            continue;
        }

        // Parse the data
        if (sscanf(line, "%lf %lf %lf %lf",
            &wavelengths[rowNum][0], &wavelengths[rowNum][1], 
            &wavelengths[rowNum][2], &wavelengths[rowNum][3]) == 4)
        {
            // Go to the next row
            rowNum++;
        }

    }

    fclose(fp);

    return 0;
}

double eV_2_micron(double energy_in_eV)
{

    /*
    
    This function calculates the wavelength in micrometers that corresponds to the energy inputted.

    */

    double h = 4.1357e-15;   // eV s
    double c = 3e8;          // m/s
    
    double wavelength = (h * c / energy_in_eV) * 1e6;    // in microns

    return wavelength;
}


int get_wavelengths(double lambda_min, double lambda_max, double **wavelengths, int numRows_wavelengths, double *center_wavelengths, int *wavelengths_indices)
{
    /*
    This function returns the center wavelengths, wavelength indices and the size of these arrays. Center wavelengths stores the Habing radiation field wavelengths in micrometer, wavelength indices
    stores the indices of the emissivity values for radiationaField_J.dat file corresponding to Habing convention wavelengths. The wavelength_indices array stores not the indices of the 
    wavelength values inside the radiationField_wavelengths.dat file but it stores the indices on the emissivity values that corresponds to the Habing wavelengths in the radiationField_J.dat file.
    */

    int estimatedSize = numRows_wavelengths; // Assuming maximum possible size for simplicity.

    int counter = 0;
    for (int i = 0; i < numRows_wavelengths; i++) 
    {
        if (wavelengths[i][0] <= lambda_max && wavelengths[i][0] >= lambda_min) 
        {
            center_wavelengths[counter] = wavelengths[i][0];
            wavelengths_indices[counter] = i + 1;
            counter++;
        }
    }
    return counter;
}


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

double isrf_calculator_in_Habing_units(double *emissivity_data, int *wavelengths_indices, int size_of_wavelength_indices, double *center_wavelengths)
{

    // Create an array to filter the emissivity data. 
    double *y = malloc(sizeof(double) * size_of_wavelength_indices);
    if (y == NULL)
    {
        perror("Unable to allocate memory for the filtered data.");
        exit(1);
    }

    // Set up the y array such that it only contains emissivity values corresponding to the Habing radiation field wavelengths
    for (int i = 0; i < size_of_wavelength_indices; i++)
    {
        // Filter the data such that it will only contain the Habing radiation field wavelengths
        y[i] = emissivity_data[wavelengths_indices[i]];
    }

    // Calculate the Habing radiation field in Habing Units. Eqn 1.9 in Rybicki and Lightman and Eqn 12.6, Draine...
    // Integrate the angle averaged emissivity values
    double integrated_emissivity = simpson_integral(y, center_wavelengths, size_of_wavelength_indices);     // W/m^2
    // Calculating the total radiation density (J/m^3). Eqn 1.9 in Rybicki and Lightman
    double u = integrated_emissivity * 4 * pi / c;                                                          // J/m^3
    // Find G parameter
    double G = u / u_HAB;                                                                                   // Unitless
    // Normalizing the total radiation density with habing radiation field to find the isrf in Habing units. Eqn 12.6, Draine
    double G_over_G0 = G / G0;                                                                              // In Habing Units

    // Free the allocated memory
    free(y);

    return G_over_G0;
}


int free_memory(int numRows, double **data)
{
    // Free dynamically allocated memory
    for (int i = 0; i < numRows; i++) 
    {
        free(data[i]);
    }
    free(data);
}

int main() 
{
    // Start the counter: 
    double timer = 0;  

    clock_t begin = clock();

    // Code starts here 


    // Read the radiationField_J.dat file
    // Define the name of the file
    const char *filename = "FB15N1024_gal5_z0_grid_radiationField_J.dat"; 

    int numRows_radiationField = countDataRows(filename);
    
    // Dynamically allocate memory for radiationField based on numRows_radiationField
    double **radiationField = malloc(numRows_radiationField * sizeof(double*));
    for (int i = 0; i < numRows_radiationField; i++) {
        radiationField[i] = malloc(NUM_COLUMNS_J * sizeof(double));
    }

    readRadiationField_J_dat_File(filename, radiationField);


    // Read the *radiationField_wavelengths.dat file 
    const char *filename_wavelengths = "FB15N1024_gal5_z0_grid_radiationField_wavelengths.dat";

    // Get the number of rows in that file
    int numRows_wavelengths = countDataRows(filename_wavelengths);

    //Dynamically allocate memory 
    double **wavelengths = malloc(numRows_wavelengths * sizeof(double*));
    for (int i = 0; i < numRows_wavelengths; i++)
    {
        wavelengths[i] = malloc(NUM_COLUMNS_WAVELENGTHS * sizeof(double*));
    }

    // Read the file 
    readWWavelenghts(filename_wavelengths, wavelengths);

/////////////////////////////////////////////////////////////////////////////////////////////////////

    /*
    TODO: Delete below.
    This part is for printing the read data
    */

    // Print the radiationField array (for testing)
    for (int i = 0; i < 10 && i < numRows_radiationField; i++) {
        for (int j = 0; j < NUM_COLUMNS_J; j++) {
            printf("%e ", radiationField[i][j]);
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < numRows_wavelengths; i++)
    {
        for (int j = 0; j < NUM_COLUMNS_WAVELENGTHS; j++)
        {
            printf("%e ", wavelengths[i][j]);
        }
        printf("\n");
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////

    // Habing radiation field is defined between the 6 eV and 13.6 eV
    double lambda_min = eV_2_micron(13.6); 
    double lambda_max = eV_2_micron(6);

    // Allocate memory to store the center_wavelengths and the indices of the these wavelengths in the radiationField_J.dat file.
    double *center_wavelengths = (double *)malloc(numRows_wavelengths * sizeof(double));
    int *wavelengths_indices = (int *)malloc(numRows_wavelengths * sizeof(int));
    if (!center_wavelengths || !wavelengths_indices) {
        // Handle memory allocation error
        perror("Unable to allocate memory!");
        exit(1);
    }

    int size_of_wavelength_indices = get_wavelengths(lambda_min, lambda_max, wavelengths, numRows_wavelengths, center_wavelengths, wavelengths_indices);

    // Trim the arrays to the actual size if needed
    // Be careful! I saw that if the new size is smaller than the old size then realloc might result with the loss of data as it works by truncating the data. My new allocation is indeed smaller than 
    // the original one, but I have checked and the values are correct. 
    center_wavelengths = (double *)realloc(center_wavelengths, size_of_wavelength_indices * sizeof(double));
    wavelengths_indices = (int *) realloc(wavelengths_indices, size_of_wavelength_indices * sizeof(int));    

    // Calculate the Habing radiation field in the Habing units.
    double *G_over_G0 = (double *)malloc(numRows_radiationField * sizeof(double));
    if (G_over_G0 == NULL)
    {
        // Handel memory allocation error
        perror("Unable to allocate memory!");
        exit(1);
    }

    for (int i = 0; i < numRows_radiationField; i++)
    {
        G_over_G0[i] = isrf_calculator_in_Habing_units(radiationField[i], wavelengths_indices, size_of_wavelength_indices, center_wavelengths);
    }

    

/////////////////////////////////////////////////////////////////////////////////////////////////////

    // Free the allocated memories.
    free_memory(numRows_radiationField, radiationField);
    free_memory(numRows_wavelengths, wavelengths);
    free(center_wavelengths);
    free(wavelengths_indices);
    free(G_over_G0);

    clock_t end = clock();
    // calculate elapsed time by finding difference (end - begin) and
    // dividing the difference by CLOCKS_PER_SEC to convert to seconds
    timer += (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Time it took for data to run is: %f\n", timer);

    return 0;
}
