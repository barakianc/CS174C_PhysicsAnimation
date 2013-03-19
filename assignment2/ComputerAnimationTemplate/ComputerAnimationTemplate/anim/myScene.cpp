////////////////////////////////////////////////////
// // Template code for  CS 174C
////////////////////////////////////////////////////

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
#include <string.h>
#include <util/util.h>
#include <GLModel/GLModel.h>
#include "anim.h"
#include "animTcl.h"
#include "myScene.h"
#include "SampleParticle.h"
#include "SampleGravitySimulator.h"
#include <util/jama/tnt_stopwatch.h>
#include <util/jama/jama_lu.h>

// register a sample variable with the shell.
// Available types are:
// - TCL_LINK_INT 
// - TCL_LINK_FLOAT

int g_testVariable = 10;

SETVAR myScriptVariables[] = {
	"testVariable", TCL_LINK_INT, (char *) &g_testVariable,
	"",0,(char *) NULL
};


//---------------------------------------------------------------------------------
//			Hooks that are called at appropriate places within anim.cpp
//---------------------------------------------------------------------------------

//this is the particle class it holds all the information about a particle
struct particle{
	VectorObj p_pos;
	VectorObj p_velocity;
	VectorObj p_force;
	VectorObj p_prevpos;
	double p_mass;
	bool p_lock;
};
//spring class holds the parameters for a spring
struct spring{
	int s_indx1,s_indx2;
	double s_ks;
	double s_kd;
	double s_restlngth;
};
//this is the system that contains the particles
class PartSys : public BaseSystem{
public:
	PartSys(const std::string& name);
	~PartSys();
	virtual void getState(particle *p);
	virtual void setState(particle *p);
	void display(GLenum mode = GL_RENDER);
	int command(int argc,char **argv);
	int getNumPart();
	int getNumSpr();
	particle* getParts();
	spring* getSprings();
	void initSprings(int snum);
	bool setSpring(int in1,int in2, double s_ks,double s_kd,double s_rlngth);
private:
	particle *p_list;
	spring *sim_sprlist;
	int num_parts;
	int num_springs;
	int num_setsprngs;
};
//the simulator that controls the particles
class PartSim : public BaseSimulator{
public:
	PartSim(const std::string& name, PartSys* target);
	~PartSim();
	int step(double time);
	int init(double time){
		return 0;
	};
	int command(int argc,char **argv);
	void euler(particle *p);
	void symplectic(particle *p);
	void verlet(particle *p);
	void computeNewForces();
	void pressedSpring(int p_x, int p_y);
	void released();
private:
	PartSys* prtsys;
	double sym_ground_ks;
	double sym_ground_kd;
	double sym_grav;
	double sym_kdrag;
	double sim_time_step;
	int sim_num_springs;
	int sim_integ;
	bool sim_verlet_frst;
	bool sim_button_pressed;
	int click_closest;
	int clk_x,clk_y;
};

//gloabl PartSim pointer allows for the mymouse function to tell the particle simulator
//when the left button is pressed down and released
PartSim* g_sim1;
// start or end interaction
void myMouse(int button, int state, int x, int y)
{

	// let the global resource manager know about the new state of the mouse 
	// button
	GlobalResourceManager::use()->setMouseButtonInfo( button, state );

	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
	{
		//call the pressedSpring() function when the left button is pressed
		Vector res;
		pickFromXYPlane(res,x,y);
		g_sim1->pressedSpring(res[0],res[1]);
		/*
		animTcl::OutputMessage(
			"My mouse received a mouse button press event\n");
		*/
	}
	if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
	{
		//tell the particle simulator when the mouse button is released
		g_sim1->released();
		/*
		animTcl::OutputMessage(
			"My mouse received a mouse button release event\n") ;
		*/
	}
}	// myMouse

// interaction (mouse motion)
void myMotion(int x, int y)
{

	GLMouseButtonInfo updatedMouseButtonInfo = 
		GlobalResourceManager::use()->getMouseButtonInfo();

	if( updatedMouseButtonInfo.button == GLUT_LEFT_BUTTON )
	{	/*
		animTcl::OutputMessage(
			"My mouse motion callback received a mousemotion event\n") ;
			*/
	}

}	// myMotion


void MakeScene(void)
{
	bool success;

	//make partSys
	PartSys* part1 = new PartSys("partSys");
	success = GlobalResourceManager::use()->addSystem(part1,true);
	assert(success);

	//make partSim
	g_sim1 = new PartSim("partSim",NULL);
	success = GlobalResourceManager::use()->addSimulator(g_sim1);
	assert(success);

	/* 
	
	This is where you instantiate all objects, systems, and simulators and 
	register them with the global resource manager
	*/

	/* SAMPLE SCENE */
	/*
	bool success;

	// register a system
	SampleParticle* sphere1 = new SampleParticle( "sphere1" );

	success = GlobalResourceManager::use()->addSystem( sphere1, true );

	// make sure it was registered successfully
	assert( success );

	// register a simulator
	SampleGravitySimulator* gravitySimulator = 
		new SampleGravitySimulator( "gravity", sphere1 );

	success = GlobalResourceManager::use()->addSimulator( gravitySimulator );

	// make sure it was registered successfully
	assert( success );
	*/
	/* END SAMPLE SCENE */

	// the following code shows you how to retrieve a system that was registered 
	// with the resource manager. 
	/*
	BaseSystem* sampleSystemRetrieval;

	// retrieve the system
	sampleSystemRetrieval = 
		GlobalResourceManager::use()->getSystem( "sphere1" );

	// make sure you got it
	assert( sampleSystemRetrieval );

	BaseSimulator* sampleSimulatorRetrieval;

	// retrieve the simulator
	sampleSimulatorRetrieval = 
		GlobalResourceManager::use()->getSimulator( "gravity" );

	// make sure you got it
	assert( sampleSimulatorRetrieval );
	*/
}	// MakeScene

// OpenGL initialization
void myOpenGLInit(void)
{
	animTcl::OutputMessage("Initialization routine was called.");

}	// myOpenGLInit

void myIdleCB(void)
{
	
	return;

}	// myIdleCB

void myKey(unsigned char key, int x, int y)
{
	 animTcl::OutputMessage("My key callback received a key press event\n");
	return;

}	// myKey

int testGlobalCommand(ClientData clientData, Tcl_Interp *interp, int argc,	char **argv)
{
	 animTcl::OutputMessage("This is a test command!");
	return TCL_OK;

}	// testGlobalCommand

void mySetScriptCommands(Tcl_Interp *interp)
{

	// here you can register additional generic (they do not belong to any object) 
	// commands with the shell

	Tcl_CreateCommand(interp, "test", testGlobalCommand, (ClientData) NULL,
					  (Tcl_CmdDeleteProc *)	NULL);

}	// mySetScriptCommands
//PartSys constructor
PartSys::PartSys(const std::string &name): BaseSystem( name){
	this->num_parts = 0;
	this->num_springs = 0;
	this->p_list = NULL;
	this->sim_sprlist = NULL;
}
//PartSys destructor
PartSys::~PartSys(){
	delete [] this->p_list;
	delete [] this->sim_sprlist;
}
//impliments the paticle system's commands
int PartSys::command(int argc, char **argv){
	if( argc < 1 )
	{
		animTcl::OutputMessage("system %s: wrong number of params.", m_name) ;
		return TCL_ERROR ;
	}
	//impliments the dim command
	else if( strcmp(argv[0],"dim")==0){
		if(argc == 2){
			this->num_parts = (int)atof(argv[1]);
			this->p_list = new particle[this->num_parts];
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of inputs!");
			return TCL_ERROR;
		}
	}
	//impliments the particle command
	else if(strcmp(argv[0], "particle")==0){
		if(argc == 9){
			int indx = (int)atof(argv[1]);
			this->p_list[indx].p_mass = atof(argv[2]);
			this->p_list[indx].p_pos[0] = atof(argv[3]);
			this->p_list[indx].p_pos[1] = atof(argv[4]);
			this->p_list[indx].p_pos[2] = atof(argv[5]);
			this->p_list[indx].p_velocity[0] = atof(argv[6]);
			this->p_list[indx].p_velocity[1] = atof(argv[7]);
			this->p_list[indx].p_velocity[2] = atof(argv[8]);
			this->p_list[indx].p_lock = false;
			this->p_list[indx].p_force[0] = 0.0;
			this->p_list[indx].p_force[1] = 0.0;
			this->p_list[indx].p_force[2] = 0.0;
			this->p_list[indx].p_prevpos[0] = 0.0;
			this->p_list[indx].p_prevpos[1] = 0.0;
			this->p_list[indx].p_prevpos[2] = 0.0;
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of parameters!");
			return TCL_ERROR;
		}
	}
	//impliments the all_velocities command
	else if(strcmp(argv[0],"all_velocities") == 0){
		if(argc == 4){
			for(int i =0; i < this->num_parts; i++){
				this->p_list[i].p_velocity[0] = atof(argv[1]);
				this->p_list[i].p_velocity[1] = atof(argv[2]);
				this->p_list[i].p_velocity[2] = atof(argv[3]);
			}
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of parameters!");
			return TCL_ERROR;
		}
	}
}
//getState and setState are not used
void PartSys::getState(particle *p){
}
void PartSys::setState(particle *p){
}
//draw the plane, particles, and springs
void PartSys::display(GLenum mode){
	glPushMatrix();
	glScalef(.125,.125,.125);
	glEnable(GL_LIGHTING);
	//glEnable(GL_TEXTURE_2D);
	glDisable (GL_BLEND); 
	glDisable (GL_DITHER); 
	glDisable (GL_FOG); 
	glDisable (GL_LIGHTING); 
	glDisable (GL_TEXTURE_1D); 
	glDisable (GL_TEXTURE_2D); 
	//glDisable (GL_TEXTURE_3D); 
	glShadeModel (GL_FLAT);
	glMatrixMode(GL_MODELVIEW);

	//draw the xz plane
	glPushMatrix() ;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glBegin(GL_QUADS);
	glColor3f(0.5f,0.5f,0.0f);
	glVertex3f(50.0f,0.0f,50.0f);
	glColor3f(0.5f,0.5f,0.0f);
	glVertex3f(-50.0f,0.0f,50.0f);
	glColor3f(0.5f,0.5f,0.0f);
	glVertex3f(-50.0f,0.0f,-50.0f);
	glColor3f(0.5f,0.5f,0.0f);
	glVertex3f(50.0f,0.0f,-50.0f);
	glPopMatrix();
	glPopAttrib();
	glEnd();

	//draw all the particles

	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPointSize(6.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0f,1.0f,1.0f);

	for(int i =0; i < this->num_parts; i++){
		Vector res;
		//pickFromXYPlane(res,(int)(),(int)());res[0]res[1]
		glVertex3f((float)(this->p_list[i].p_pos[0]),(float)(this->p_list[i].p_pos[1]),(float)(this->p_list[i].p_pos[2]));
	}
	glEnd();
	glPopMatrix();
	glPopAttrib();

	//draw all Springs
	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glLineWidth(2.0f);
	glBegin(GL_LINES);
	for(int i = 0; i<this->num_springs; i++){
		glColor3f(1.0f,0.0f,0.0f);
		Vector res;
		//pickFromXYPlane(res,(int)(),(int)());res[0]res[1]
		glVertex3f((float)(this->p_list[this->sim_sprlist[i].s_indx1].p_pos[0]),(float)(this->p_list[this->sim_sprlist[i].s_indx1].p_pos[1]),(float)(this->p_list[this->sim_sprlist[i].s_indx1].p_pos[2]));
		glColor3f(1.0f,0.0f,0.0f);
		//pickFromXYPlane(res,(int)(),(int)());res[0]res[1]
		glVertex3f((float)(this->p_list[this->sim_sprlist[i].s_indx2].p_pos[0]),(float)(this->p_list[this->sim_sprlist[i].s_indx2].p_pos[1]),(float)(this->p_list[this->sim_sprlist[i].s_indx2].p_pos[2]));
	}
	glEnd();
	glPopMatrix();
	glPopAttrib();

	glPopMatrix();
}
//get the number of particles on the system
int PartSys::getNumPart(){
	return this->num_parts;
}
//get the number of springs in the system
int PartSys::getNumSpr(){
	return this->num_springs;
}
//initialize all springs
void PartSys::initSprings(int snum){
	this->sim_sprlist = new spring[snum];
	this->num_springs = snum;
	this->num_setsprngs = 0;
	for(int i=0; i<snum; i++){
		this->sim_sprlist[i].s_indx1=0;
		this->sim_sprlist[i].s_indx2=0;
		this->sim_sprlist[i].s_kd = 0.0;
		this->sim_sprlist[i].s_ks = 0.0;
		this->sim_sprlist[i].s_restlngth = 0.0;
	}
}
//sets a spring to the given input values
bool PartSys::setSpring(int in1,int in2,double s_ks,double s_kd,double s_rlngth){
	if( this->num_setsprngs < this->num_springs){
		this->sim_sprlist[this->num_setsprngs].s_indx1 = in1;
		this->sim_sprlist[this->num_setsprngs].s_indx2 = in2;
		this->sim_sprlist[this->num_setsprngs].s_ks = s_ks;
		this->sim_sprlist[this->num_setsprngs].s_kd = s_kd;
		this->sim_sprlist[this->num_setsprngs].s_restlngth = s_rlngth;
		this->num_setsprngs++;
		return true;
	}
	else{
		return false;
	}
}
particle* PartSys::getParts(){
	return this->p_list;
}
spring* PartSys::getSprings(){
	return this->sim_sprlist;
}
//PartSim constructor
PartSim::PartSim(const std::string &name, PartSys *target):BaseSimulator( name )
{

	this->prtsys = target;
	this->sim_time_step = .001;
	this->sim_verlet_frst = false;
	this->sym_grav = -9.81;
	this->sym_ground_ks = 1000.0;
	this->sym_ground_kd = 100.0;
	this->sym_kdrag = 0.0;
	this->sim_integ = 0;
	this->sim_button_pressed = false;
}
PartSim::~PartSim(){
}
//impliments the step class for the particle simulator
int PartSim::step(double time){
	//first compute the force on each of the particels
	this->computeNewForces();

	particle *partlst = this->prtsys->getParts();
	int numparts = this->prtsys->getNumPart();
	//now for each particle use the specified integration method
	for(int i = 0; i < numparts; i++){
		if(this->sim_integ == 0){
			this->euler(&(partlst[i]));
		}
		else if(this->sim_integ == 1){
			this->symplectic(&(partlst[i]));
		}
		else if(this->sim_integ == 2){
			this->verlet(&(partlst[i]));
		}
	}
	if(this->sim_integ == 2){
		this->sim_verlet_frst = false;
	}
	return 0;
}
int PartSim::command(int argc, char **argv){
	if( argc < 1 )
	{
		animTcl::OutputMessage("system %s: wrong number of params.", m_name) ;
		return TCL_ERROR ;
	}
	//impliment the link command
	else if(strcmp(argv[0],"link") == 0){
		if(argc == 3){
			BaseSystem* parts;

			// retrieve the system
			parts = GlobalResourceManager::use()->getSystem( argv[1]);

			// make sure you got it
			assert( parts );

			this->prtsys = (PartSys*)parts;
			this->sim_num_springs = (int)atof(argv[2]);
			this->prtsys->initSprings(this->sim_num_springs);
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of input parameters");
			return TCL_ERROR;
		}
	}
	//impliment the spring command
	else if(strcmp(argv[0],"spring") == 0){
		if(argc == 6){
			int indx1 = (int)atof(argv[1]);
			int indx2 = (int)atof(argv[2]);
			double spr_ks = atof(argv[3]);
			double spr_kd = atof(argv[4]);
			double spr_rlngth = atof(argv[5]);
			bool success = this->prtsys->setSpring(indx1,indx2,spr_ks,spr_kd,spr_rlngth);
			if(success){
				return TCL_OK;
			}
			else{
				animTcl::OutputMessage("Could not set spring, reached specified limit of springs");
				return TCL_OK;
			}
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of input elements!");
			return TCL_ERROR;
		}
	}
	//impliment the fix command
	else if(strcmp(argv[0],"fix") == 0){
		if(argc == 2){
			int pindx = (int)atof(argv[1]);
			particle *part = this->prtsys->getParts();
			part[pindx].p_lock = true;
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("Error: Incorrect number of inputs!");
			return TCL_ERROR;
		}
	}
	//impliment the integration command
	else if(strcmp(argv[0],"integration") == 0){
		if(argc == 3){
			if(strcmp(argv[1],"euler")==0){
				this->sim_integ = 0;
				this->sim_time_step = atof(argv[2]);
			}
			else if(strcmp(argv[1],"symplectic")==0){
				this->sim_integ = 1;
				this->sim_time_step = atof(argv[2]);
			}
			else if(strcmp(argv[1],"verlet")==0){
				this->sim_integ = 2;
				this->sim_verlet_frst = true;
				this->sim_time_step = atof(argv[2]);
			}
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of inputs!");
			return TCL_ERROR;
		}
	}
	//impliment the ground command
	else if(strcmp(argv[0],"ground")==0){
		if(argc == 3){
			this->sym_ground_ks = atof(argv[1]);
			this->sym_ground_kd = atof(argv[2]);
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of inputs!");
			return TCL_ERROR;
		}
	}
	//impliment the gravity command
	else if(strcmp(argv[0],"gravity") == 0){
		if(argc == 2){
			this->sym_grav = atof(argv[1]);
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of Inputs!");
			return TCL_ERROR;
		}
	}
	//impliment the grag command
	else if(strcmp(argv[0],"drag") == 0){
		if(argc == 2){
			this->sym_kdrag = atof(argv[1]);
			return TCL_OK;
		}
		else{
			animTcl::OutputMessage("ERROR: Incorrect number of inputs!");
			return TCL_ERROR;
		}
	}

}
void PartSim::computeNewForces(){
	int numparts = this->prtsys->getNumPart();
	particle *partlst = this->prtsys->getParts();
	spring *sprnglst = this->prtsys->getSprings();
	//first zero out all previous forces
	for(int i=0; i< numparts; i++){
		partlst[i].p_force[0] = 0.0;
		partlst[i].p_force[1] = 0.0;
		partlst[i].p_force[2] = 0.0;
	}
	//compute force of gravity for each particle
	for(int i =0; i< numparts; i++){
		partlst[i].p_force[1] += ((this->sym_grav)*(partlst[i].p_mass));
	}

	//compute force from clicked point
	if(this->sim_button_pressed){
		Vector clked;
		clked[0] = this->clk_x;
		clked[1] = this->clk_y;
		clked[2] = partlst[this->click_closest].p_pos[2];
		VectorObj dir = (partlst[this->click_closest].p_pos) - clked;
		double mag = dir.length();
		VectorObj nrml = dir.normalize();
		VectorObj veloc = partlst[this->click_closest].p_velocity;
		VectorObj cl_fsp = nrml * ((1000)*(-1)*(mag));
		VectorObj cl_fd = nrml * ((-1)*(10)*(veloc.dot(nrml)));
		partlst[this->click_closest].p_force += cl_fsp;
		partlst[this->click_closest].p_force += cl_fd;
	}
	//compute dampering force
	for(int i = 0; i< numparts; i++){
		partlst[i].p_force[0] += ((-1*(this->sym_kdrag))*(partlst[i].p_velocity[0]));
		partlst[i].p_force[1] += ((-1*(this->sym_kdrag))*(partlst[i].p_velocity[1]));
		partlst[i].p_force[2] += ((-1*(this->sym_kdrag))*(partlst[i].p_velocity[2]));
	}
	//compute the force on each particle for each spring
	for(int i = 0; i <this->sim_num_springs; i++){
		if(sprnglst[i].s_indx1 != sprnglst[i].s_indx2){
			VectorObj vectr = (partlst[sprnglst[i].s_indx2].p_pos)-(partlst[sprnglst[i].s_indx1].p_pos);
			double veclngth = vectr.length();
			VectorObj nrmlvect = vectr.normalize();
			VectorObj veloc = (partlst[sprnglst[i].s_indx2].p_velocity) - (partlst[sprnglst[i].s_indx1].p_velocity);
			VectorObj fsp = nrmlvect*((sprnglst[i].s_ks)*((sprnglst[i].s_restlngth)-(veclngth)));
			VectorObj fd = nrmlvect*((-1)*(sprnglst[i].s_kd)*((veloc.dot(nrmlvect))));
			partlst[sprnglst[i].s_indx2].p_force += fsp;
			partlst[sprnglst[i].s_indx2].p_force += fd;
			partlst[sprnglst[i].s_indx1].p_force += (-1)*fsp;
			partlst[sprnglst[i].s_indx1].p_force += (-1)*fd;
		}
	}

	//apply the force to particles below the y axis
	for(int i = 0; i< numparts; i++){
		VectorObj normal;
		normal[0] = 0.0;
		normal[1] = 1.0;
		normal[2] = 0.0;
		if(partlst[i].p_pos[1] < 0.0){
			VectorObj gfsp = (normal)*((this->sym_ground_ks)*(partlst[i].p_pos.dot(normal)));
			VectorObj gfd = normal* ((-1)*(this->sym_ground_kd)*(partlst[i].p_velocity.dot(normal)));
			partlst[i].p_force -= gfsp;
			partlst[i].p_force += gfd;
		}
	}

}
//impliments euler forward integration
void PartSim::euler(particle *p){
	if(!(p->p_lock)){
		p->p_prevpos[0] = p->p_pos[0];
		p->p_prevpos[1] = p->p_pos[1];
		p->p_prevpos[2] = p->p_pos[2];
		p->p_pos[0] = p->p_pos[0] + ((this->sim_time_step)*(p->p_velocity[0]));
		p->p_pos[1] = p->p_pos[1] + ((this->sim_time_step)*(p->p_velocity[1]));
		p->p_pos[2] = p->p_pos[2] + ((this->sim_time_step)*(p->p_velocity[2]));
		p->p_velocity[0] = p->p_velocity[0] + ((this->sim_time_step)*((p->p_force[0])/(p->p_mass)));
		p->p_velocity[1] = p->p_velocity[1] + ((this->sim_time_step)*((p->p_force[1])/(p->p_mass)));
		p->p_velocity[2] = p->p_velocity[2] + ((this->sim_time_step)*((p->p_force[2])/(p->p_mass)));
	}
}
//impliments symplectic euler integration
void PartSim::symplectic(particle *p){
	if(!(p->p_lock)){
		p->p_prevpos[0] = p->p_pos[0];
		p->p_prevpos[1] = p->p_pos[1];
		p->p_prevpos[2] = p->p_pos[2];
		p->p_velocity[0] = p->p_velocity[0] + ((this->sim_time_step)*((p->p_force[0])/(p->p_mass)));
		p->p_velocity[1] = p->p_velocity[1] + ((this->sim_time_step)*((p->p_force[1])/(p->p_mass)));
		p->p_velocity[2] = p->p_velocity[2] + ((this->sim_time_step)*((p->p_force[2])/(p->p_mass)));
		p->p_pos[0] = p->p_pos[0] + ((this->sim_time_step)*(p->p_velocity[0]));
		p->p_pos[1] = p->p_pos[1] + ((this->sim_time_step)*(p->p_velocity[1]));
		p->p_pos[2] = p->p_pos[2] + ((this->sim_time_step)*(p->p_velocity[2]));
	}
}
//impliments verlet integration
void PartSim::verlet(particle *p){
	if(!(p->p_lock)){
		if(this->sim_verlet_frst){
			//this->sim_verlet_frst = false;
			this->symplectic(p);
		}
		else{
			double prevposx = p->p_pos[0];
			double prevposy = p->p_pos[1];
			double prevposz = p->p_pos[2];
			p->p_pos[0] = (2*(p->p_pos[0]))-(p->p_prevpos[0]) + (((this->sim_time_step)*(this->sim_time_step))*((p->p_force[0])/(p->p_mass)));
			p->p_pos[1] = (2*(p->p_pos[1]))-(p->p_prevpos[1]) + (((this->sim_time_step)*(this->sim_time_step))*((p->p_force[1])/(p->p_mass)));
			p->p_pos[2] = (2*(p->p_pos[2]))-(p->p_prevpos[2]) + (((this->sim_time_step)*(this->sim_time_step))*((p->p_force[2])/(p->p_mass)));
			p->p_velocity[0] = p->p_velocity[0] + ((this->sim_time_step)*((p->p_force[0])/(p->p_mass)));
			p->p_velocity[1] = p->p_velocity[1] + ((this->sim_time_step)*((p->p_force[1])/(p->p_mass)));
			p->p_velocity[2] = p->p_velocity[2] + ((this->sim_time_step)*((p->p_force[2])/(p->p_mass)));
			p->p_prevpos[0] = prevposx;
			p->p_prevpos[1] = prevposy;
			p->p_prevpos[2] = prevposz;
		}
	}

}

void PartSim::released(){
	this->sim_button_pressed = false;
}
//handles when the user clicks the screen and determines the closest particle
void PartSim::pressedSpring(int p_x, int p_y){
	this->sim_button_pressed = true;
	particle *p = this->prtsys->getParts();
	this->clk_x = p_x+10;
	this->clk_y = p_y+10;
	VectorObj clicked;
	clicked[0] = p_x+10;
	clicked[1] = p_y+10;
	double minlngth;
	int temp_min;
	for(int i =0; i< (this->prtsys->getNumPart()); i++){
		clicked[2] = p[i].p_pos[2];
		VectorObj lngth = clicked - p[i].p_pos;
		if(i == 0){
			temp_min = i;
			minlngth = lngth.length();
		}
		else{
			if( (lngth.length() < minlngth)){
				temp_min = i;
				minlngth = lngth.length();
			}
		}
	}
	this->click_closest = temp_min;

}