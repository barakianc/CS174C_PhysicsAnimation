#ifdef WIN32
#include <windows.h>
#endif


#include <GL/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef WIN32
#include <GL/glut.h>
#endif

#ifdef WIN32
#include "GL/glut.h"
#endif

#include <shared/defs.h>

#include <util/util.h>
#include "anim.h"
#include <tcl/tcl.h>
#include "animTcl.h"

#include "myScene.h"

#include "GlobalResourceManager.h"

void DrawScene(GLenum mode) 
{
	//for( int i = 0 ; i < NumObjects ; i++ )
	//	Objects[i]->display(mode) ;

	GlobalResourceManager::use()->display( mode );

}

void InitSimulation(void)
{
	// call the init of all simulators
		GlobalResourceManager::use()->initializeAllSimulators();

}

void SimulationStep(void)
{

	GlobalResourceManager::use()->advanceSimulationTime();

	GlobalResourceManager::use()->stepAllSimulators();

	// animTcl::OutputMessage("This is the simulation time %lf", SIMTIME);
}
