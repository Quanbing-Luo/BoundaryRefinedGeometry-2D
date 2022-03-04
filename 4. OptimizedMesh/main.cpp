#include <iostream>
#include <string>
#include "OptimizedMesh.h"


int main()
{		
	OptimizedMesh m;
	m.input();	
	m.initialization();
	m.optimization();
	m.output();
}



