#include "InitialMesh.h"

int main()
{		
	
	InitialMesh im;

	im.input();	

	im.initialization();
	
	im.meshGeneration(4.0, 6.0, 1.0, 1.1);

	im.output();	
}



