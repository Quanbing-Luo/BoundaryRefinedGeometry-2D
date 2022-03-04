#include <iostream>
#include <string>
#include "Mesh.h"

int main()
{		
	
	Mesh m;

	m.input();

	m.initialization();	

	m.meshGeneration();

	m.output();	
}



