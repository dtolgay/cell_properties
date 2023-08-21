#include <stdio.h>

int main()
{
	double dummy = 6.266159e-09;

	double changed; 

	changed = dummy * 1e12; 

	printf("%f\n", changed);

	return 0;
}