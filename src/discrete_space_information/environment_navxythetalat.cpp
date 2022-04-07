/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include <cstring>
#include <ctime>
#include <sbpl/discrete_space_information/environment_navxythetalat.h>
#include <sbpl/utils/2Dgridsearch.h>
#include <sbpl/utils/key.h>
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>

#if TIME_DEBUG
static clock_t time3_addallout = 0;
static clock_t time_gethash = 0;
static clock_t time_createhash = 0;
static clock_t time_getsuccs = 0;
#endif

static long int checks = 0;

#define XYTHETA2INDEX(X,Y,THETA) (THETA + X*EnvNAVXYTHETALATCfg.NumThetaDirs + \
                                  Y*EnvNAVXYTHETALATCfg.EnvWidth_c*EnvNAVXYTHETALATCfg.NumThetaDirs)

#define XYTHETACART2INDEX(X,Y,THETA,CARTANGLE) (THETA + X*EnvNAVXYTHETALATCartCfg.NumThetaDirs + \
                                  Y*EnvNAVXYTHETALATCartCfg.EnvWidth_c*EnvNAVXYTHETALATCartCfg.NumThetaDirs + CARTANGLE*EnvNAVXYTHETALATCartCfg.EnvWidth_c*EnvNAVXYTHETALATCartCfg.EnvHeight_c*EnvNAVXYTHETALATCartCfg.NumThetaDirs)


EnvironmentNAVXYTHETALATTICE::EnvironmentNAVXYTHETALATTICE()
{
    grid2Dsearchfromstart = NULL;
    grid2Dsearchfromgoal = NULL;
    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;
    iteration = 0;
    bucketsize = 0; // fixed bucket size
    blocksize = 1;
    bUseNonUniformAngles = false;

    // for trailer EnvCfg
    EnvNAVXYTHETALATCartCfg.obsthresh = ENVNAVXYTHETALATCART_DEFAULTOBSTHRESH;
	EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh = EnvNAVXYTHETALATCartCfg.obsthresh; //the value that pretty much makes it disabled
	EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh = -1; //the value that pretty much makes it disabled

	EnvNAVXYTHETALATCART.bInitialized = false;

	EnvNAVXYTHETALATCartCfg.actionwidth = NAVXYTHETALATCART_DEFAULT_ACTIONWIDTH;

    EnvNAVXYTHETALATCartCfg.NumThetaDirs = NAVXYTHETALATCART_THETADIRS;


	//no memory allocated in cfg yet
	EnvNAVXYTHETALATCartCfg.Grid2D = NULL;
	EnvNAVXYTHETALATCartCfg.ActionsV = NULL;
	EnvNAVXYTHETALATCartCfg.PredActionsV = NULL;
}

EnvironmentNAVXYTHETALATTICE::~EnvironmentNAVXYTHETALATTICE()
{
    SBPL_PRINTF("destroying XYTHETALATTICE\n");
    if (grid2Dsearchfromstart != NULL) {
        delete grid2Dsearchfromstart;
    }
    grid2Dsearchfromstart = NULL;

    if (grid2Dsearchfromgoal != NULL) {
        delete grid2Dsearchfromgoal;
    }
    grid2Dsearchfromgoal = NULL;

    // for trailer EnvCfg
    if(EnvNAVXYTHETALATCartCfg.Grid2D != NULL)
	{	
		for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) 
			delete [] EnvNAVXYTHETALATCartCfg.Grid2D[x];
		delete [] EnvNAVXYTHETALATCartCfg.Grid2D;
		EnvNAVXYTHETALATCartCfg.Grid2D = NULL;
	}

	//delete actions
	if(EnvNAVXYTHETALATCartCfg.ActionsV != NULL)
	{
		for(int tind = 0; tind < EnvNAVXYTHETALATCartCfg.NumThetaDirs; tind++){
            for (int cind = 0; cind < CART_THETADIRS; cind++) {
                delete [] EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind];
            }
            delete [] EnvNAVXYTHETALATCartCfg.ActionsV[tind];
        }
		delete [] EnvNAVXYTHETALATCartCfg.ActionsV;
		EnvNAVXYTHETALATCartCfg.ActionsV = NULL;
	}
	if(EnvNAVXYTHETALATCartCfg.PredActionsV != NULL)
	{
		delete [] EnvNAVXYTHETALATCartCfg.PredActionsV;
		EnvNAVXYTHETALATCartCfg.PredActionsV = NULL;
	}
}

static unsigned int inthash(unsigned int key)
{
    key += (key << 12);
    key ^= (key >> 22);
    key += (key << 4);
    key ^= (key >> 9);
    key += (key << 10);
    key ^= (key >> 2);
    key += (key << 7);
    key ^= (key >> 12);
    return key;
}

void EnvironmentNAVXYTHETALATTICE::SetConfiguration(
    int width, int height, 
    const unsigned char* mapdata,
    int startx, int starty, int starttheta, int startcartangle,
    int goalx, int goaly, int goaltheta, int goalcartangle,
    double cellsize_m,
    double nominalvel_mpersecs,
    double timetoturn45degsinplace_secs,
    const std::vector<sbpl_2Dpt_t>& robot_perimeterV,
    const std::vector<sbpl_2Dpt_t> & cart_perimeterV,
    const sbpl_2Dpt_t &cart_offset)
    // const sbpl_2Dpt_t &cart_cp_offset,
    // const double rod_length)
{
    EnvNAVXYTHETALATCartCfg.EnvWidth_c = width;
    EnvNAVXYTHETALATCartCfg.EnvHeight_c = height;
    EnvNAVXYTHETALATCartCfg.StartX_c = startx;
    EnvNAVXYTHETALATCartCfg.StartY_c = starty;
    EnvNAVXYTHETALATCartCfg.StartTheta = starttheta;
    EnvNAVXYTHETALATCartCfg.StartCartAngle = startcartangle;

    if (EnvNAVXYTHETALATCartCfg.StartX_c < 0 ||
        EnvNAVXYTHETALATCartCfg.StartX_c >= EnvNAVXYTHETALATCartCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.StartY_c < 0 ||
        EnvNAVXYTHETALATCartCfg.StartY_c >= EnvNAVXYTHETALATCartCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.StartTheta < 0 ||
        EnvNAVXYTHETALATCartCfg.StartTheta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates for theta");
    }
    if (EnvNAVXYTHETALATCartCfg.StartCartAngle < 0 || 
       EnvNAVXYTHETALATCartCfg.StartCartAngle >= CART_THETADIRS) 
    {
        throw SBPL_Exception("ERROR: illegal start coordinates for theta cart");
    }

    EnvNAVXYTHETALATCartCfg.EndX_c = goalx;
    EnvNAVXYTHETALATCartCfg.EndY_c = goaly;
    EnvNAVXYTHETALATCartCfg.EndTheta = goaltheta;
    EnvNAVXYTHETALATCartCfg.EndCartAngle = goalcartangle;

    if (EnvNAVXYTHETALATCartCfg.EndX_c < 0 ||
        EnvNAVXYTHETALATCartCfg.EndX_c >= EnvNAVXYTHETALATCartCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.EndY_c < 0 ||
        EnvNAVXYTHETALATCartCfg.EndY_c >= EnvNAVXYTHETALATCartCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.EndTheta < 0 ||
        EnvNAVXYTHETALATCartCfg.EndTheta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates for theta");
    }
    if (EnvNAVXYTHETALATCartCfg.EndCartAngle < 0 || 
        EnvNAVXYTHETALATCartCfg.EndCartAngle >= CART_THETADIRS) 
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates for theta cart");
    }

    EnvNAVXYTHETALATCartCfg.FootprintPolygon = robot_perimeterV;
    EnvNAVXYTHETALATCartCfg.CartPolygon = cart_perimeterV;
    EnvNAVXYTHETALATCartCfg.CartOffset = cart_offset;

    EnvNAVXYTHETALATCartCfg.CartCenterOffset.x = 0.0;
    EnvNAVXYTHETALATCartCfg.CartCenterOffset.y = 0.0;

    for(unsigned int i= 0; i < cart_perimeterV.size(); i++)
    {
        EnvNAVXYTHETALATCartCfg.CartCenterOffset.x += cart_perimeterV[i].x;
        EnvNAVXYTHETALATCartCfg.CartCenterOffset.y += cart_perimeterV[i].y;
    }

    EnvNAVXYTHETALATCartCfg.CartCenterOffset.x = EnvNAVXYTHETALATCartCfg.CartCenterOffset.x/cart_perimeterV.size();
    EnvNAVXYTHETALATCartCfg.CartCenterOffset.y = EnvNAVXYTHETALATCartCfg.CartCenterOffset.y/cart_perimeterV.size();


    EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs = nominalvel_mpersecs;
    EnvNAVXYTHETALATCartCfg.cellsize_m = cellsize_m;
    EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs = timetoturn45degsinplace_secs;

    // unallocate the 2D environment
    if (EnvNAVXYTHETALATCartCfg.Grid2D != NULL) {
        for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
            delete[] EnvNAVXYTHETALATCartCfg.Grid2D[x];
        }
        delete[] EnvNAVXYTHETALATCartCfg.Grid2D;
        EnvNAVXYTHETALATCartCfg.Grid2D = NULL;
    }

    // allocate the 2D environment
    EnvNAVXYTHETALATCartCfg.Grid2D = new unsigned char*[EnvNAVXYTHETALATCartCfg.EnvWidth_c];
    for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
        EnvNAVXYTHETALATCartCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETALATCartCfg.EnvHeight_c];
    }

    // environment:
    if (0 == mapdata) {
        for (int y = 0; y < EnvNAVXYTHETALATCartCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
                EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = 0;
            }
        }
    }
    else {
        for (int y = 0; y < EnvNAVXYTHETALATCartCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
                EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = mapdata[x + y * width];
            }
        }
    }
}

void EnvironmentNAVXYTHETALATTICE::ReadConfiguration(FILE* fCfg)
{
    // read in the configuration of environment and initialize
    // EnvNAVXYTHETALATCfg structure
    char sTemp[1024], sTemp1[1024];
    int dTemp;
    int x, y;

    // discretization(cells)
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    strcpy(sTemp1, "discretization(cells):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format (discretization)" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    EnvNAVXYTHETALATCartCfg.EnvWidth_c = atoi(sTemp);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    EnvNAVXYTHETALATCartCfg.EnvHeight_c = atoi(sTemp);

    // Scan for optional NumThetaDirs parameter. Check for following obsthresh.
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "NumThetaDirs:");
    if (strcmp(sTemp1, sTemp) != 0) {
        // optional NumThetaDirs not available; default is NAVXYTHETALAT_THETADIRS (16)
        strcpy(sTemp1, "obsthresh:");
        if (strcmp(sTemp1, sTemp) != 0) {
            std::stringstream ss;
            ss << "ERROR: configuration file has incorrect format" <<
                    " Expected " << sTemp1 << " got " << sTemp;
            throw SBPL_Exception(ss.str());
        }
        else {
            EnvNAVXYTHETALATCartCfg.NumThetaDirs = NAVXYTHETALAT_THETADIRS;
        }
    }
    else {
        if (fscanf(fCfg, "%s", sTemp) != 1) {
            throw SBPL_Exception("ERROR: ran out of env file early (NumThetaDirs)");
        }
        EnvNAVXYTHETALATCartCfg.NumThetaDirs = atoi(sTemp);

        //obsthresh:
        if (fscanf(fCfg, "%s", sTemp) != 1) {
            throw SBPL_Exception("ERROR: ran out of env file early (obsthresh)");
        }
        strcpy(sTemp1, "obsthresh:");
        if (strcmp(sTemp1, sTemp) != 0) {
            std::stringstream ss;
            ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
            throw SBPL_Exception(ss.str());
        }
    }

    // obsthresh
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.obsthresh = atoi(sTemp);
    SBPL_PRINTF("obsthresh = %d\n", EnvNAVXYTHETALATCartCfg.obsthresh);

    //cost_inscribed_thresh:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cost_inscribed_thresh:");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh = atoi(sTemp);
    SBPL_PRINTF("cost_inscribed_thresh = %d\n", EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh);

    //cost_possibly_circumscribed_thresh:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cost_possibly_circumscribed_thresh:");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh = atoi(sTemp);
    SBPL_PRINTF("cost_possibly_circumscribed_thresh = %d\n", EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh);

    //cellsize
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cellsize(meters):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.cellsize_m = atof(sTemp);

    //speeds
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "nominalvel(mpersecs):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
            " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs = atof(sTemp);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "timetoturn45degsinplace(secs):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs = atof(sTemp);

    // start(meters,rads):
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.StartX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETALATCartCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.StartY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETALATCartCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }

    EnvNAVXYTHETALATCartCfg.StartTheta_rad = atof(sTemp);

    if (EnvNAVXYTHETALATCartCfg.StartX_c < 0 ||
        EnvNAVXYTHETALATCartCfg.StartX_c >= EnvNAVXYTHETALATCartCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.StartY_c < 0 ||
        EnvNAVXYTHETALATCartCfg.StartY_c >= EnvNAVXYTHETALATCartCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }

    // end(meters,rads):
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.EndX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETALATCartCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvNAVXYTHETALATCartCfg.EndY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETALATCartCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }

    EnvNAVXYTHETALATCartCfg.EndTheta_rad = atof(sTemp);

    if (EnvNAVXYTHETALATCartCfg.EndX_c < 0 ||
        EnvNAVXYTHETALATCartCfg.EndX_c >= EnvNAVXYTHETALATCartCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal end coordinates");
    }
    if (EnvNAVXYTHETALATCartCfg.EndY_c < 0 ||
        EnvNAVXYTHETALATCartCfg.EndY_c >= EnvNAVXYTHETALATCartCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal end coordinates");
    }

    // unallocate the 2d environment
    if (EnvNAVXYTHETALATCartCfg.Grid2D != NULL) {
        for (x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
            delete[] EnvNAVXYTHETALATCartCfg.Grid2D[x];
        }
        delete[] EnvNAVXYTHETALATCartCfg.Grid2D;
        EnvNAVXYTHETALATCartCfg.Grid2D = NULL;
    }

    // allocate the 2D environment
    EnvNAVXYTHETALATCartCfg.Grid2D = new unsigned char*[EnvNAVXYTHETALATCartCfg.EnvWidth_c];
    for (x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
        EnvNAVXYTHETALATCartCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETALATCartCfg.EnvHeight_c];
    }

    // environment:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    for (y = 0; y < EnvNAVXYTHETALATCartCfg.EnvHeight_c; y++) {
        for (x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
            if (fscanf(fCfg, "%d", &dTemp) != 1) {
                throw SBPL_Exception("ERROR: incorrect format of config file");
            }
            EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = dTemp;
        }
    }
}

// trailer
bool EnvironmentNAVXYTHETALATTICE::ReadinCell(
    EnvNAVXYTHETALATCART3Dcell_t* cell,
    FILE* fIn)
{
   char sTemp[60];

	if(fscanf(fIn, "%s", sTemp) == 0){
	   return false;
    }
	cell->x = atoi(sTemp);
	if(fscanf(fIn, "%s", sTemp) == 0){
	   return false;
    }
	cell->y = atoi(sTemp);
	if(fscanf(fIn, "%s", sTemp) == 0){
	   return false;
    }
	cell->theta = atoi(sTemp);
    //normalize the angle TODO Yucheng: necessary? why only theta not cart theta?
	cell->theta = NORMALIZEDISCTHETA(cell->theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);

	if(fscanf(fIn, "%s", sTemp) == 0)
	   return false;
	cell->cartangle = atoi(sTemp);
    //normalize the angle
    // TODO 16.03
    // cell-
  //	cell->cartangle = NORMALIZEDISCTHETA(cell->cartangle, CART_THETADIRS);

	return true;
}

bool EnvironmentNAVXYTHETALATTICE::ReadinPose(
    EnvNAVXYTHETALATCART3Dpt_t* pose,
    FILE* fIn)
{
    char sTemp[60];

    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->x = atof(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->y = atof(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->theta = atof(sTemp);
    pose->theta = normalizeAngle(pose->theta);

    if(fscanf(fIn, "%s", sTemp) == 0){
	   return false;
    }
	pose->cartangle = atof(sTemp);

    return true;
}

int EnvironmentNAVXYTHETALATTICE::normalizeDiscAngle(int theta) const
{
    if (bUseNonUniformAngles) {
        if (theta < 0) {
            theta += EnvNAVXYTHETALATCartCfg.NumThetaDirs;
        }
        if (theta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs) {
            theta -= EnvNAVXYTHETALATCartCfg.NumThetaDirs;
        }
    }
    else {
        theta = NORMALIZEDISCTHETA(theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
    }
    return theta;
}

double EnvironmentNAVXYTHETALATTICE::DiscTheta2ContNew(int theta) const
{
    if (bUseNonUniformAngles) {
        return DiscTheta2ContFromSet(theta);
    }
    else {
        return DiscTheta2Cont(theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
    }
}

int EnvironmentNAVXYTHETALATTICE::ContTheta2DiscNew(double theta) const
{
    if (bUseNonUniformAngles) {
        return ContTheta2DiscFromSet(theta);
    }
    else {
        return ContTheta2Disc(theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
    }
}

// trailer angle conversion
double CartDiscTheta2Cont(int nTheta, int NUMOFANGLEVALS)
{
  double thetaBinSize = (2*MAX_CART_ANGLE)/(CART_THETADIRS-1);
  return -MAX_CART_ANGLE+nTheta*thetaBinSize;
//   double thetaBinSize = (2*180/(double)NUMOFANGLEVALS);
//   // SBPL_INFO("What is it? %d, %f, %d",nTheta, thetaBinSize, NUMOFANGLEVALS);
//   double result = (double)nTheta*thetaBinSize;
//   // SBPL_INFO("What is it? %f", result);

//   return result;
}

//converts continuous (radians) version of angle into discrete
//maps 0->0, [delta/2, 3/2*delta)->1, [3/2*delta, 5/2*delta)->2,...
int CartContTheta2Disc(double fTheta, int NUMOFANGLEVALS)
{
    // for 360° version
    // SBPL_INFO("cont phi to discrete value %f, %d", (fTheta+(2*M_PI/NUMOFANGLEVALS)/2)/(2*M_PI/NUMOFANGLEVALS), floor((fTheta+(2*M_PI/NUMOFANGLEVALS)/2)/(2*M_PI/NUMOFANGLEVALS)));
    // return floor((fTheta+(2*M_PI/NUMOFANGLEVALS)/2)/(2*M_PI/NUMOFANGLEVALS));
  
  // for +-45° version
  if(fTheta < -MAX_CART_ANGLE)
    return 0;
  else if (fTheta > MAX_CART_ANGLE)
    return (CART_THETADIRS-1);
  else
  {
    double thetaBinSize = (2*MAX_CART_ANGLE)/(CART_THETADIRS-1);
    // SBPL_INFO("cont phi to discrete value %f, %d",fTheta, (int)((fTheta+MAX_CART_ANGLE+thetaBinSize/2.0)/thetaBinSize));
    
    return (int)((fTheta+MAX_CART_ANGLE+thetaBinSize/2.0)/thetaBinSize);
  }
}

double EnvironmentNAVXYTHETALATTICE::DiscTheta2ContFromSet(int theta) const
{
    theta = normalizeDiscAngle(theta);

    // ThetaDirs should contain extra angle (2PI) for overlap
    if (EnvNAVXYTHETALATCartCfg.NumThetaDirs >= (int)EnvNAVXYTHETALATCartCfg.ThetaDirs.size()) {
        throw SBPL_Exception("ERROR: list of bin angles are not properly set to use function DiscTheta2ConfFromSet");
    }

    if (theta > EnvNAVXYTHETALATCartCfg.NumThetaDirs || theta < 0) {
        std::stringstream ss;
        ss << "ERROR: discrete value theta " << theta << " out of range";
        throw SBPL_Exception(ss.str());
    }
    return EnvNAVXYTHETALATCartCfg.ThetaDirs[theta];
}

int EnvironmentNAVXYTHETALATTICE::ContTheta2DiscFromSet(double theta) const
{
    theta = normalizeAngle(theta);
    // ThetaDirs should contain extra angle (2PI) for overlap
    if (EnvNAVXYTHETALATCartCfg.NumThetaDirs >= (int) EnvNAVXYTHETALATCartCfg.ThetaDirs.size()) {
        throw SBPL_Exception("ERROR: list of bin angles are not properly set to use function ContTheta2DiscFromSet");
    }

    int lower_bound_ind = -1;
    int upper_bound_ind = -1;
    for (int i = 1; i < (int) EnvNAVXYTHETALATCartCfg.ThetaDirs.size(); i++) {
        if ((EnvNAVXYTHETALATCartCfg.ThetaDirs[i]) >= theta) {
            lower_bound_ind = i - 1;
            upper_bound_ind = i;
            break;
        }
    }

    // Critical error if could not find bin location from given angle
    if (lower_bound_ind == -1) {
        std::stringstream ss;
        ss << "ERROR: unable to find bin index for angle " << theta;
        throw SBPL_Exception(ss.str());
    }

    // Get closest angle of two
    double angle_low = EnvNAVXYTHETALATCartCfg.ThetaDirs[lower_bound_ind];
    double angle_up = EnvNAVXYTHETALATCartCfg.ThetaDirs[upper_bound_ind];
    double diff_low = fabs(theta - angle_low);
    double diff_up = fabs(theta - angle_up);

    if (diff_low < diff_up) {
        return lower_bound_ind;
    }
    else {
        // Wrap upper bound index around when it reaches last index (assumed to be 2PI)
        if (upper_bound_ind == EnvNAVXYTHETALATCartCfg.NumThetaDirs) {
            upper_bound_ind = 0;
        }
        return upper_bound_ind;
    }
}

// reconstruct function for trailer
bool EnvironmentNAVXYTHETALATTICE::ReadinMotionPrimitive(
    SBPL_xythetacart_mprimitive* pMotPrim,
    FILE* fIn)
{
    char sTemp[1024];
    int dTemp;
    char sExpected[1024];
    int numofIntermPoses;
    float fTemp;

    // read in actionID
    strcpy(sExpected, "primID:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        fflush(stdout);
        return false;
    }
    if (fscanf(fIn, "%d", &pMotPrim->motprimID) != 1) {
        return false;
    }

    // read in start angle
    strcpy(sExpected, "startangle_c:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &dTemp) == 0) {
        SBPL_ERROR("ERROR reading startangle\n");
        return false;
    }
    pMotPrim->starttheta_c = dTemp;

    // read in start trailer angle
    strcpy(sExpected, "startjointangle_c:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &dTemp) == 0) {
        SBPL_ERROR("ERROR reading startjointangle\n");
        return false;
    }
    pMotPrim->startphi_c = dTemp;

    // read in end pose
    strcpy(sExpected, "endpose_c:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }

    if (ReadinCell(&pMotPrim->endcell, fIn) == false) {
        SBPL_ERROR("ERROR: failed to read in endsearchpose\n");
        return false;
    }

    // read in action cost
    strcpy(sExpected, "additionalactioncostmult:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &dTemp) != 1) {
        return false;
    }
    pMotPrim->additionalactioncostmult = dTemp;

    if (bUseNonUniformAngles) {
        // read in action turning radius
        strcpy(sExpected, "turning_radius:");
        if (fscanf(fIn, "%s", sTemp) == 0) {
            return false;
        }
        if (strcmp(sTemp, sExpected) != 0) {
            SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
            return false;
        }
        if (fscanf(fIn, "%f", &fTemp) != 1) {
            return false;
        }
        pMotPrim->turning_radius = fTemp;
    }

    // read in intermediate poses
    strcpy(sExpected, "intermediateposes:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &numofIntermPoses) != 1) {
        return false;
    }
    // all intermposes should be with respect to 0,0 as starting pose since it
    // will be added later and should be done after the action is rotated by
    // initial orientation
    for (int i = 0; i < numofIntermPoses; i++) {
        EnvNAVXYTHETALATCART3Dpt_t intermpose;
        if (ReadinPose(&intermpose, fIn) == false) {
            SBPL_ERROR("ERROR: failed to read in intermediate poses\n");
            return false;
        }
        pMotPrim->intermptV.push_back(intermpose);
    }

    // Check that the last pose of the motion matches (within lattice
    // resolution) the designated end pose of the primitive
    EnvNAVXYTHETALATCART3Dpt_t sourcepose;
    sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
    sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
    sourcepose.theta = DiscTheta2ContNew(pMotPrim->starttheta_c);
    // SBPL_INFO("startphi_c!!!!!%d", pMotPrim->startphi_c);
    sourcepose.cartangle = CartDiscTheta2Cont(pMotPrim->startphi_c, CART_THETADIRS);
    // SBPL_INFO("cartangle!!!!!%f, %f", sourcepose.cartangle, CartDiscTheta2Cont(pMotPrim->startphi_c, CART_THETADIRS));



    double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
    double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
    double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;
    double mp_endcartangle_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].cartangle;

    int endx_c = CONTXY2DISC(mp_endx_m, EnvNAVXYTHETALATCartCfg.cellsize_m);
    int endy_c = CONTXY2DISC(mp_endy_m, EnvNAVXYTHETALATCartCfg.cellsize_m);
    int endtheta_c = ContTheta2DiscNew(mp_endtheta_rad);
    int endcartangle_c = CartContTheta2Disc(mp_endcartangle_rad, CART_THETADIRS);
    // SBPL_INFO("end pose %d, %d, %d, %d. movement primitive pose %d, %d, %d, %d\n", 
    //           endx_c, endy_c, endtheta_c, endcartangle_c, pMotPrim->endcell.x, pMotPrim->endcell.y, pMotPrim->endcell.theta, pMotPrim->endcell.cartangle);
    if (endx_c != pMotPrim->endcell.x ||
        endy_c != pMotPrim->endcell.y ||
        endtheta_c != pMotPrim->endcell.theta ||
        endcartangle_c != pMotPrim->endcell.cartangle)
    {
        SBPL_ERROR( "ERROR: incorrect primitive %d with startangle=%d and startcartangle=%d"
                   "last interm point %f %f %f %f does not match end pose %d %d %d %d\n",
                   pMotPrim->motprimID, pMotPrim->starttheta_c, pMotPrim->startphi_c,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].cartangle,
                   pMotPrim->endcell.x, 
                   pMotPrim->endcell.y,
                   pMotPrim->endcell.theta,
                   pMotPrim->endcell.cartangle);
        SBPL_FFLUSH(stdout);
        return false;
    }

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::ReadMotionPrimitivesCart(FILE* fMotPrims)
{
    char sTemp[1024], sExpected[1024];
    float fTemp;
    int dTemp;
    int totalNumofActions = 0;


    SBPL_INFO("Reading in motion primitives...");
    fflush(stdout);

    //read in the resolution
    strcpy(sExpected, "resolution_m:");
    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        fflush(stdout);
        return false;
    }
    if (fscanf(fMotPrims, "%f", &fTemp) == 0) {
        return false;
    }
    if (fabs(fTemp - EnvNAVXYTHETALATCartCfg.cellsize_m) > ERR_EPS) {
        SBPL_ERROR("ERROR: invalid resolution %f (instead of %f) in the dynamics file\n", fTemp, EnvNAVXYTHETALATCartCfg.cellsize_m);
        fflush(stdout);
        return false;
    }
    SBPL_INFO("resolution_m: %f\n", fTemp);

    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    SBPL_INFO("sTemp: %s\n", sTemp);
    if (strncmp(sTemp, "min_turning_radius_m:", 21) == 0) {
        bUseNonUniformAngles = true;
    }
    SBPL_INFO("bUseNonUniformAngles = %d", bUseNonUniformAngles);

    if (bUseNonUniformAngles) {
        float min_turn_rad;
        strcpy(sExpected, "min_turning_radius_m:");
        if (strcmp(sTemp, sExpected) != 0) {
            SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
            fflush(stdout);
            return false;
        }
        if (fscanf(fMotPrims, "%f", &min_turn_rad) == 0) {
            return false;
        }
        SBPL_PRINTF("min_turn_rad: %f\n", min_turn_rad);
        fflush(stdout);
        if (fscanf(fMotPrims, "%s", sTemp) == 0) {
            return false;
        }
    }

    // read in the angular resolution
    strcpy(sExpected, "numberofangles:");
    if (strcmp(sTemp, sExpected) != 0) {
       SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
       return false;
    }
    if (fscanf(fMotPrims, "%d", &dTemp) == 0) {
        return false;
    }
    if (dTemp != EnvNAVXYTHETALATCartCfg.NumThetaDirs) {
        SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n", dTemp, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
        return false;
    }
    SBPL_PRINTF("numberofangles: %d\n", dTemp);
    EnvNAVXYTHETALATCartCfg.NumThetaDirs = dTemp;

    if (bUseNonUniformAngles) {
        // read in angles
        EnvNAVXYTHETALATCartCfg.ThetaDirs.clear();
        for (int i = 0; i < EnvNAVXYTHETALATCartCfg.NumThetaDirs; i++)
        {
            std::ostringstream string_angle_index;
            string_angle_index << i;
            std::string angle_string = "angle:" + string_angle_index.str();

            float angle;
            strcpy(sExpected, angle_string.c_str());
            if (fscanf(fMotPrims, "%s", sTemp) == 0) {
                return false;
            }
            if (strcmp(sTemp, sExpected) != 0) {
                SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
                return false;
            }
            if (fscanf(fMotPrims, "%f", &angle) == 0) {
                return false;
            }
            SBPL_PRINTF("%s %f\n", angle_string.c_str(), angle);
            EnvNAVXYTHETALATCartCfg.ThetaDirs.push_back(angle);
        }
        EnvNAVXYTHETALATCartCfg.ThetaDirs.push_back(2.0 * M_PI); // Add 2 PI at end for overlap
    }

    // modified for trailer angle, not configured in the cpp file
    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    SBPL_INFO("sTemp: %s\n", sTemp);
    strcpy(sExpected, "numberofjointangles:");
    if (strcmp(sTemp, sExpected) != 0) {
       SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
       return false;
    }
    if (fscanf(fMotPrims, "%d", &dTemp) == 0) {
        return false;
    }
    SBPL_PRINTF("numberofjointangles: %d\n", dTemp);

    // TODO Yucheng: add NumThetaDirs for trailer!!!
    if (dTemp != CART_THETADIRS) {
        SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n", dTemp, CART_THETADIRS);
        return false;
    }

    // read in the total number of actions
    strcpy(sExpected, "totalnumberofprimitives:");
    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fMotPrims, "%d", &totalNumofActions) == 0) {
        return false;
    }
    SBPL_PRINTF("totalnumberofprimitives: %d\n", totalNumofActions);

    // Read in motion primitive for each action
    for (int i = 0; i < totalNumofActions; i++) {
        SBPL_xythetacart_mprimitive motprim;

        if (!EnvironmentNAVXYTHETALATTICE::ReadinMotionPrimitive(&motprim, fMotPrims)) {
            return false;
        }

        EnvNAVXYTHETALATCartCfg.mprimV.push_back(motprim);

    }
    SBPL_PRINTF("done");
    SBPL_FFLUSH(stdout);
    return true;
}

void EnvironmentNAVXYTHETALATTICE::ComputeReplanningDataforAction(
    EnvNAVXYTHETALATCARTAction_t* action)
{
    int j;

    // iterate over all the cells involved in the action
    EnvNAVXYTHETALATCART3Dcell_t startcell3d, endcell3d;
    for (int i = 0; i < (int)action->intersectingcellsV.size(); i++) {
        // SBPL_INFO("intersectingcellsV: %d", (int)action->intersectingcellsV.size());
        // compute the translated affected search Pose - what state has an
        // outgoing action whose intersecting cell is at 0,0
        startcell3d.theta = action->starttheta;
        startcell3d.cartangle = action->startcartangle;
        startcell3d.x = -action->intersectingcellsV.at(i).x;
        startcell3d.y = -action->intersectingcellsV.at(i).y;

        // SBPL_INFO("action: x,y (%d, %d), theta (%d, %d), cartangle (%d, %d)", action->dX, action->dY, action->starttheta, action->endtheta, action->startcartangle, action->endcartangle);
        // SBPL_INFO("startcell: (%d, %d, %d, %d)", startcell3d.x, startcell3d.y, startcell3d.theta, startcell3d.cartangle);
        // compute the translated affected search Pose - what state has an
        // incoming action whose intersecting cell is at 0,0
        if (bUseNonUniformAngles) {
            endcell3d.theta = normalizeDiscAngle(action->endtheta);
        }
        else {
            endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
        }
        endcell3d.x = startcell3d.x + action->dX;
        endcell3d.y = startcell3d.y + action->dY;
        endcell3d.cartangle = action->endcartangle;
        // SBPL_INFO("endcell: (%d, %d, %d, %d)", endcell3d.x, endcell3d.y, endcell3d.theta, endcell3d.cartangle);

        //store the cells if not already there
        for (j = 0; j < (int)affectedsuccstatesVc.size(); j++) {
            if (affectedsuccstatesVc.at(j) == endcell3d) {
                break;
            }
        }
        if (j == (int)affectedsuccstatesVc.size()) {
            affectedsuccstatesVc.push_back(endcell3d);
        }

        for (j = 0; j < (int)affectedpredstatesVc.size(); j++) {
            if (affectedpredstatesVc.at(j) == startcell3d) {
                break;
            }
        }
        if (j == (int)affectedpredstatesVc.size()) {
            affectedpredstatesVc.push_back(startcell3d);
        }
    } // over intersecting cells

    // add the centers since with h2d we are using these in cost computations
    // ---intersecting cell = origin
    // compute the translated affected search Pose - what state has an outgoing
    // action whose intersecting cell is at 0,0
    startcell3d.theta = action->starttheta;
    startcell3d.x = -0;
    startcell3d.y = -0;
    // startcell3d.cartangle = CartContTheta2Disc(0.0, CART_THETADIRS);
    startcell3d.cartangle = action->startcartangle;

    // compute the translated affected search Pose - what state has an incoming
    // action whose intersecting cell is at 0,0
    if (bUseNonUniformAngles) {
        endcell3d.theta = normalizeDiscAngle(action->endtheta);
    }
    else {
        endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
    }
    endcell3d.x = startcell3d.x + action->dX;
    endcell3d.y = startcell3d.y + action->dY;
    // endcell3d.cartangle = CartContTheta2Disc(0.0, CART_THETADIRS);
    endcell3d.cartangle =NORMALIZEDISCTHETA(action->endcartangle, CART_THETADIRS);


    //store the cells if not already there
    for (j = 0; j < (int)affectedsuccstatesVc.size(); j++) {
        if (affectedsuccstatesVc.at(j) == endcell3d) {
            break;
        }
    }
    if (j == (int)affectedsuccstatesVc.size()) {
        affectedsuccstatesVc.push_back(endcell3d);
    }

    for (j = 0; j < (int)affectedpredstatesVc.size(); j++) {
        if (affectedpredstatesVc.at(j) == startcell3d) {
            break;
        }
    }
    if (j == (int)affectedpredstatesVc.size()) {
        affectedpredstatesVc.push_back(startcell3d);
    }

    //---intersecting cell = outcome state
    // compute the translated affected search Pose - what state has an outgoing
    // action whose intersecting cell is at 0,0
    startcell3d.theta = action->starttheta;
    startcell3d.x = -action->dX;
    startcell3d.y = -action->dY;
    // 22.03 added
    startcell3d.cartangle = action->startcartangle;

    // compute the translated affected search Pose - what state has an incoming
    // action whose intersecting cell is at 0,0
    if (bUseNonUniformAngles) {
        endcell3d.theta = normalizeDiscAngle(action->endtheta);
    }
    else {
        endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
    }
    endcell3d.x = startcell3d.x + action->dX;
    endcell3d.y = startcell3d.y + action->dY;
    // added
    endcell3d.cartangle =NORMALIZEDISCTHETA(action->endcartangle, CART_THETADIRS);

    for (j = 0; j < (int)affectedsuccstatesVc.size(); j++) {
        if (affectedsuccstatesVc.at(j) == endcell3d) {
            break;
        }
    }
    if (j == (int)affectedsuccstatesVc.size()) {
        affectedsuccstatesVc.push_back(endcell3d);
    }

    for (j = 0; j < (int)affectedpredstatesVc.size(); j++) {
        if (affectedpredstatesVc.at(j) == startcell3d) {
            break;
        }
    }
    if (j == (int)affectedpredstatesVc.size()) {
        affectedpredstatesVc.push_back(startcell3d);
    }
}

void EnvironmentNAVXYTHETALATTICE::ComputeReplanningDataCart()
{
    // iterate over all actions
    // orientations
    for (int tind = 0; tind < EnvNAVXYTHETALATCartCfg.NumThetaDirs; tind++) {
        // cartangles
        for (int cind = 0; cind < CART_THETADIRS; cind++) {
            // actions
            for (int aind = 0; aind < EnvNAVXYTHETALATCartCfg.actionwidth; aind++) {
                // compute replanning data for this action
                SBPL_INFO("compute replanning Data for actions! tind %d, cind %d, aind %d\n", tind,cind, aind);
                ComputeReplanningDataforAction(&EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind]);
            }
        }
    }
}

// for trailer
void EnvironmentNAVXYTHETALATTICE::PrecomputeActionswithCompleteMotionPrimitive(
    std::vector<SBPL_xythetacart_mprimitive>* motionprimitiveV)
{
    SBPL_PRINTF("Pre-computing action data using motion primitives for every angle...\n");
    // EnvNAVXYTHETALATCartCfg.ActionsV = new EnvNAVXYTHETALATCARTAction_t*[EnvNAVXYTHETALATCartCfg.NumThetaDirs];
    
    EnvNAVXYTHETALATCartCfg.ActionsV = new EnvNAVXYTHETALATCARTAction_t**[EnvNAVXYTHETALATCartCfg.NumThetaDirs];

    EnvNAVXYTHETALATCartCfg.PredActionsV = new std::vector<EnvNAVXYTHETALATCARTAction_t*>[EnvNAVXYTHETALATCartCfg.NumThetaDirs];
    std::vector<sbpl_2Dcell_t> footprint;

    if (motionprimitiveV->size() % EnvNAVXYTHETALATCartCfg.NumThetaDirs != 0) {
        throw SBPL_Exception("ERROR: motionprimitives should be uniform across actions");
    }

    EnvNAVXYTHETALATCartCfg.actionwidth = ((int)motionprimitiveV->size()) / (EnvNAVXYTHETALATCartCfg.NumThetaDirs*CART_THETADIRS);

    // iterate over source angles
    int maxnumofactions = 0;
    for (int tind = 0; tind < EnvNAVXYTHETALATCartCfg.NumThetaDirs; tind++) {
        SBPL_PRINTF("pre-computing for angle %d out of %d angles\n", tind, EnvNAVXYTHETALATCartCfg.NumThetaDirs);
        
        // EnvNAVXYTHETALATCartCfg.ActionsV[tind] = new EnvNAVXYTHETALATCARTAction_t[EnvNAVXYTHETALATCartCfg.actionwidth];
        EnvNAVXYTHETALATCartCfg.ActionsV[tind] = new EnvNAVXYTHETALATCARTAction_t*[CART_THETADIRS];

        //Yucheng: iterate over source cartangle (phi)
        for (int cind = 0; cind < CART_THETADIRS; cind++)
        // change to the member variable of EnvNAVXYTHETALATCartCfg.Num
        {
            EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind] = new EnvNAVXYTHETALATCARTAction_t[EnvNAVXYTHETALATCartCfg.actionwidth];
        
            // compute sourcepose
            EnvNAVXYTHETALATCART3Dpt_t sourcepose;
            sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
            sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
            sourcepose.theta = DiscTheta2ContNew(tind);
            sourcepose.cartangle = CartDiscTheta2Cont(cind, CART_THETADIRS);

            int numofactions = 0;
            int aind = -1;
            for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
                // SBPL_INFO("number of motion primitive: size %d,tind %d, cind %d, mind %d, aind %d", (int)motionprimitiveV->size(),tind, cind, mind, aind);

                // test, result: start phi angle didnt initialized!! 17.03
                // for (int i = 0; i < (int)motionprimitiveV->size(); i++) 
                //     SBPL_INFO("cartangle comparision, (%d, %d)", motionprimitiveV->at(i).startphi_c, cind);
                // SBPL_INFO("cartangle comparision, (%d, %d) %d", motionprimitiveV->at(mind).startphi_c, cind, motionprimitiveV->at(mind).startphi_c != cind);

                //find a motion primitive for this angle
                if (motionprimitiveV->at(mind).starttheta_c != tind || motionprimitiveV->at(mind).startphi_c != cind) {
                    continue;
                }

                aind++;
                numofactions++;

                // action index
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].aind = aind;

                // start angle
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].starttheta = tind;
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].startcartangle = cind;
        //         double mp_endx_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].x;
        // 		double mp_endy_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].y;
        // 		double mp_endtheta_rad = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].theta;
        //   //			double mp_endcartangle_rad = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].cartangle;
                
        // 		double endx = sourcepose.x + (mp_endx_m*cos(sourcepose.theta) - mp_endy_m*sin(sourcepose.theta));
        // 		double endy = sourcepose.y + (mp_endx_m*sin(sourcepose.theta) + mp_endy_m*cos(sourcepose.theta));
                
        // 		int endx_c = CONTXY2DISC(endx, EnvNAVXYTHETACARTLATCfg.cellsize_m);
        // 		int endy_c = CONTXY2DISC(endy, EnvNAVXYTHETACARTLATCfg.cellsize_m);


                // compute dislocation
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endtheta = motionprimitiveV->at(mind).endcell.theta;
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endcartangle = motionprimitiveV->at(mind).endcell.cartangle;
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dX = motionprimitiveV->at(mind).endcell.x;
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dY = motionprimitiveV->at(mind).endcell.y;

                //TODO Yucheng 17.03: modify the cost function for a reasonable calculation
                //compute cost
                if(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dY != 0 || EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dX != 0)
                    EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].cost = (int)(ceil(NAVXYTHETALAT_COSTMULT_MTOMM*EnvNAVXYTHETALATCartCfg.cellsize_m/EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs*
                                    sqrt((double)(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dX*EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dX + 
                                    EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dY*EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dY))));
                else //cost of turn in place
                    EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].cost = (int)(NAVXYTHETALAT_COSTMULT_MTOMM*
                            EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs*
                            fabs(computeMinUnsignedAngleDiff(DiscTheta2Cont(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endtheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs),
                                                            DiscTheta2Cont(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].starttheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs)))/(PI_CONST/4.0));
                //use any additional cost multiplier
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;
                //

                // compute and store interm points as well as intersecting cells
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intersectingcellsV.clear();
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.clear();
                // for test
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.reserve(10000);
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].interm3DcellsV.clear();

                EnvNAVXYTHETALATCART3Dcell_t previnterm3Dcell;
                previnterm3Dcell.x = 0;
                previnterm3Dcell.y = 0;
                previnterm3Dcell.theta = 0; 
                previnterm3Dcell.cartangle = CartContTheta2Disc(0, CART_THETADIRS); // 2 for 5 number of jointangles
                SBPL_DEBUG("Motion primitive has %d intermediate points for (%d,%d)",(int)motionprimitiveV->at(mind).intermptV.size(),tind,aind);
                // SBPL_DEBUG("Cartangle is %d", previnterm3Dcell.cartangle);


                // Compute all the intersected cells for this action (intermptV and interm3DcellsV)
                for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
                    EnvNAVXYTHETALATCART3Dpt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
                    // SBPL_PRINTF("intermpt point: %f, %f, %f, %f", intermpt.x, intermpt.y, intermpt.theta, intermpt.cartangle);

                    // SBPL_INFO("maximum aind is %d, the size of intermptV %d, %d, %d, %d",EnvNAVXYTHETALATCartCfg.actionwidth/CART_THETADIRS, EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.size(), tind, cind, aind);
                    //store it (they are with reference to 0,0,stattheta (not sourcepose.x,sourcepose.y,starttheta (that is, half-bin))
                    EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.push_back(intermpt);
                    // for (unsigned int i = 1; i < EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size(); i++) {
                    //     SBPL_DEBUG("point: (%f, %f, %f, %f)\n", EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].x, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].y, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].theta, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].cartangle);
                    // }


                    // also compute the intermediate discrete cells if not there already
                    EnvNAVXYTHETALATCART3Dpt_t pose;
                    // pose.x = intermpt.x + sourcepose.x;
                    // pose.y = intermpt.y + sourcepose.y;
                    // pose.theta = intermpt.theta;
                    pose = intermpt;
                    // SBPL_DEBUG("Pose for source point is: %f %f %f %f",intermpt.x,intermpt.y,intermpt.theta,intermpt.cartangle);

                    pose.x += sourcepose.x;
                    pose.y += sourcepose.y;
                    // SBPL_DEBUG("Pose for intermediate point %d is: %f %f %f %f",pind,pose.x,pose.y,pose.theta,pose.cartangle);
                    
                    CalculateFootprintForPose(
                        EnvNAVXYTHETALATCartCfg.FootprintPolygon,
                        EnvNAVXYTHETALATCartCfg.CartPolygon,
                        &EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intersectingcellsV,
                        pose,
                        EnvNAVXYTHETALATCartCfg.cellsize_m);
                    EnvNAVXYTHETALATCART3Dpt_t cart_center_pose = getCartCenter(pose,  EnvNAVXYTHETALATCartCfg.CartCenterOffset);

                    //now also store the intermediate discretized cell if not there already
                    EnvNAVXYTHETALATCART3Dcell_t intermediate3dCell;
                    intermediate3dCell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETALATCartCfg.cellsize_m);
                    intermediate3dCell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETALATCartCfg.cellsize_m);
                    intermediate3dCell.theta = ContTheta2Disc(pose.theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs); 
                    intermediate3dCell.cartangle = CartContTheta2Disc(pose.cartangle, CART_THETADIRS); 

                    EnvNAVXYTHETALATCART3Dcell_t intermediate3dCellCart;
                    intermediate3dCellCart.x = CONTXY2DISC(cart_center_pose.x, EnvNAVXYTHETALATCartCfg.cellsize_m);
                    intermediate3dCellCart.y = CONTXY2DISC(cart_center_pose.y, EnvNAVXYTHETALATCartCfg.cellsize_m);
                    intermediate3dCellCart.theta = ContTheta2Disc(cart_center_pose.theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs); 
                    intermediate3dCellCart.cartangle = CartContTheta2Disc(cart_center_pose.cartangle, CART_THETADIRS); 

                    // add unique cells to the list
                    if (EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].interm3DcellsV.size() == 0 ||
                        previnterm3Dcell.theta != intermediate3dCell.theta || previnterm3Dcell.x != intermediate3dCell.x || previnterm3Dcell.y != intermediate3dCell.y || previnterm3Dcell.cartangle != intermediate3dCell.cartangle)
                    {
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].interm3DcellsV.push_back(intermediate3dCell);
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].interm3DcellsV.push_back(intermediate3dCellCart);
                    }

                    previnterm3Dcell = intermediate3dCell;
                    // SBPL_DEBUG("Intersecting cells have size: %d",(int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.size());
                }

                // TODO:31.03 add remove footprint
                // RemoveSourceFootprint(sourcepose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);

                SBPL_DEBUG("action tind=%2d cind=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) endphi=%3d (%6.2f degs -> %6.2f degs) "
                        "(mprimID %3d: %3d %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
                        tind,
                        cind,
                        aind,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dX,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].dY,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endtheta,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[0].theta * 180 / PI_CONST,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.size() - 1].theta * 180 / PI_CONST,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endcartangle,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[0].cartangle * 180 / PI_CONST,
                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.size() - 1].cartangle * 180 / PI_CONST,
                        motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
                        motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta, motionprimitiveV->at(mind).endcell.cartangle,
                        (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].interm3DcellsV.size(),
                        (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intersectingcellsV.size());

                // compute linear and angular time
                // TODO Yucheng: if necessary to calculate the cartangle?
                // SBPL_DEBUG("Calculate linear and angular time now!");
                double linear_distance = 0;
                for (unsigned int i = 1; i < EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV.size(); i++) {
                    // SBPL_DEBUG("point: (%f, %f, %f, %f)\n", EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].x, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].y, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].theta, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].cartangle);
                    double x0 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[i - 1].x;
                    double y0 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[i - 1].y;
                    double x1 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[i].x;
                    double y1 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intermptV[i].y;
                    double dx = x1 - x0;
                    double dy = y1 - y0;
                    // SBPL_DEBUG("dx, dy are %f, %f. point: (%f, %f), (%f, %f)\n", dx, dy, x0, y0, x1, y1);
                    linear_distance += sqrt(dx * dx + dy * dy);
                    // SBPL_DEBUG("Linear distance is %f\n", linear_distance);
                }
                double linear_time = linear_distance / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs;
                // SBPL_DEBUG("Nominal linear velocity is %f\n", EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs);
                // SBPL_DEBUG("Linear time is %f\n", linear_time);
                double angular_distance;
                angular_distance = fabs(computeMinUnsignedAngleDiff(
                        DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endtheta),
                        DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].starttheta)));

                double angular_time = angular_distance / ((PI_CONST / 4.0) /
                                    EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs);
                // SBPL_DEBUG("Angular distance is %f, Angular time is %f\n", angular_distance, angular_time);

                // make the cost the max of the two times
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].cost =
                        (int)(ceil(NAVXYTHETALAT_COSTMULT_MTOMM * std::max(linear_time, angular_time)));
                // use any additional cost multiplier
                EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;

                // now compute the intersecting cells for this motion (including ignoring the source footprint)
                Get2dMotionCells(
                        EnvNAVXYTHETALATCartCfg.FootprintPolygon,
                        EnvNAVXYTHETALATCartCfg.CartPolygon,
                        motionprimitiveV->at(mind).intermptV,
                        &EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].intersectingcellsV,
                        EnvNAVXYTHETALATCartCfg.cellsize_m);

    // #if DEBUG
    //             SBPL_FPRINTF(fDeb,
    //                          "action tind=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) "
    //                          "cost=%4d (mprimID %3d: %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
    //                          tind,
    //                          aind,
    //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX,
    //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY,
    //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta,
    //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[0].theta * 180 / PI_CONST,
    //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size() - 1].theta * 180 / PI_CONST, EnvNAVXYTHETALATCfg.ActionsV[tind][aind].cost,
    //                          motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
    //                          motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta,
    //                          (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.size(),
    //                          (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.size());
    // #endif

                // add to the list of backward actions
                int targettheta = EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind].endtheta;
                if (targettheta < 0) {
                    targettheta = targettheta + EnvNAVXYTHETALATCartCfg.NumThetaDirs;
                }
                EnvNAVXYTHETALATCartCfg.PredActionsV[targettheta].push_back(&(EnvNAVXYTHETALATCartCfg.ActionsV[tind][cind][aind]));
            }


            if (maxnumofactions < numofactions) {
                maxnumofactions = numofactions;
            }
        }


        // compute sourcepose
        // EnvNAVXYTHETALATCART3Dpt_t sourcepose;
        // sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
        // sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETALATCartCfg.cellsize_m);
        // sourcepose.theta = DiscTheta2ContNew(tind);
        // sourcepose.cartangle = CartDiscTheta2Cont(2, CART_THETADIRS);

        // iterate over motion primitives
//         int numofactions = 0;
//         int aind = -1;
//         for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
//             SBPL_INFO("number of motion primitive: size %d, mind %d, aind %d", (int)motionprimitiveV->size(), mind, aind);

//             //find a motion primitive for this angle
//             if (motionprimitiveV->at(mind).starttheta_c != tind) {
//                 continue;
//             }

//             aind++;
//             numofactions++;

//             // action index
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].aind = aind;

//             // start angle
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].starttheta = tind;
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].startcartangle = CartContTheta2Disc(0, CART_THETADIRS);
//     //         double mp_endx_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].x;
// 	// 		double mp_endy_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].y;
// 	// 		double mp_endtheta_rad = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].theta;
//     //   //			double mp_endcartangle_rad = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size()-1].cartangle;
			
// 	// 		double endx = sourcepose.x + (mp_endx_m*cos(sourcepose.theta) - mp_endy_m*sin(sourcepose.theta));
// 	// 		double endy = sourcepose.y + (mp_endx_m*sin(sourcepose.theta) + mp_endy_m*cos(sourcepose.theta));
			
// 	// 		int endx_c = CONTXY2DISC(endx, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 	// 		int endy_c = CONTXY2DISC(endy, EnvNAVXYTHETACARTLATCfg.cellsize_m);


//             // compute dislocation
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta = motionprimitiveV->at(mind).endcell.theta;
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endcartangle = CartContTheta2Disc(0, CART_THETADIRS);
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX = motionprimitiveV->at(mind).endcell.x;
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY = motionprimitiveV->at(mind).endcell.y;

//             //
//             //compute cost
// 			if(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY != 0 || EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX != 0)
// 				EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].cost = (int)(ceil(NAVXYTHETALAT_COSTMULT_MTOMM*EnvNAVXYTHETALATCartCfg.cellsize_m/EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs*
// 								sqrt((double)(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX*EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX + 
// 								EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY*EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY))));
// 			else //cost of turn in place
// 				EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].cost = (int)(NAVXYTHETALAT_COSTMULT_MTOMM*
// 						EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs*
// 						fabs(computeMinUnsignedAngleDiff(DiscTheta2Cont(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs),
// 														DiscTheta2Cont(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].starttheta, EnvNAVXYTHETALATCartCfg.NumThetaDirs)))/(PI_CONST/4.0));
// 			//use any additional cost multiplier
// 			EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;
//             //

//             // compute and store interm points as well as intersecting cells
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.clear();
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.clear();
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.clear();

//             EnvNAVXYTHETALATCART3Dcell_t previnterm3Dcell;
//             previnterm3Dcell.x = 0;
//             previnterm3Dcell.y = 0;
//             previnterm3Dcell.theta = 0; 
//             previnterm3Dcell.cartangle = CartContTheta2Disc(0, CART_THETADIRS);
//             SBPL_DEBUG("Motion primitive has %d intermediate points for (%d,%d)",(int)motionprimitiveV->at(mind).intermptV.size(),tind,aind);
//             // SBPL_DEBUG("Cartangle is %d", previnterm3Dcell.cartangle);


//             // Compute all the intersected cells for this action (intermptV and interm3DcellsV)
//             for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
//                 EnvNAVXYTHETALATCART3Dpt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
//                 // SBPL_PRINTF("intermpt point: %f, %f, %f, %f", intermpt.x, intermpt.y, intermpt.theta, intermpt.cartangle);

//                 //store it (they are with reference to 0,0,stattheta (not sourcepose.x,sourcepose.y,starttheta (that is, half-bin))
//                 EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.push_back(intermpt);
//                 // for (unsigned int i = 1; i < EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size(); i++) {
//                 //     SBPL_DEBUG("point: (%f, %f, %f, %f)\n", EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].x, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].y, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].theta, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].cartangle);
//                 // }

//                 // also compute the intermediate discrete cells if not there already
//                 EnvNAVXYTHETALATCART3Dpt_t pose;
//                 // pose.x = intermpt.x + sourcepose.x;
//                 // pose.y = intermpt.y + sourcepose.y;
//                 // pose.theta = intermpt.theta;
//                 pose = intermpt;
//                 // SBPL_DEBUG("Pose for source point is: %f %f %f %f",intermpt.x,intermpt.y,intermpt.theta,intermpt.cartangle);

//                 pose.x += sourcepose.x;
//                 pose.y += sourcepose.y;
//                 // SBPL_DEBUG("Pose for intermediate point %d is: %f %f %f %f",pind,pose.x,pose.y,pose.theta,pose.cartangle);
                
//                 CalculateFootprintForPose(
//                     EnvNAVXYTHETALATCartCfg.FootprintPolygon,
//                     EnvNAVXYTHETALATCartCfg.CartPolygon,
//                     &EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV,
//                     pose,
//                     EnvNAVXYTHETALATCartCfg.cellsize_m);
//                 EnvNAVXYTHETALATCART3Dpt_t cart_center_pose = getCartCenter(pose,  EnvNAVXYTHETALATCartCfg.CartCenterOffset);

//                 //now also store the intermediate discretized cell if not there already
//                 EnvNAVXYTHETALATCART3Dcell_t intermediate3dCell;
//                 intermediate3dCell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETALATCartCfg.cellsize_m);
//                 intermediate3dCell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETALATCartCfg.cellsize_m);
//                 intermediate3dCell.theta = ContTheta2Disc(pose.theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs); 
// 				intermediate3dCell.cartangle = CartContTheta2Disc(pose.cartangle, CART_THETADIRS); 

//                 EnvNAVXYTHETALATCART3Dcell_t intermediate3dCellCart;
// 				intermediate3dCellCart.x = CONTXY2DISC(cart_center_pose.x, EnvNAVXYTHETALATCartCfg.cellsize_m);
// 				intermediate3dCellCart.y = CONTXY2DISC(cart_center_pose.y, EnvNAVXYTHETALATCartCfg.cellsize_m);
// 				intermediate3dCellCart.theta = ContTheta2Disc(cart_center_pose.theta, EnvNAVXYTHETALATCartCfg.NumThetaDirs); 
// 				intermediate3dCellCart.cartangle = CartContTheta2Disc(cart_center_pose.cartangle, CART_THETADIRS); 

//                 // add unique cells to the list
//                 if (EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.size() == 0 ||
//                     previnterm3Dcell.theta != intermediate3dCell.theta || previnterm3Dcell.x != intermediate3dCell.x || previnterm3Dcell.y != intermediate3dCell.y)
//                 {
//                     EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.push_back(intermediate3dCell);
//                     EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.push_back(intermediate3dCellCart);
//                 }

//                 previnterm3Dcell = intermediate3dCell;
//                 // SBPL_DEBUG("Intersecting cells have size: %d",(int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.size());
//             }

//             SBPL_DEBUG("action tind=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) "
//                        "(mprimID %3d: %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
//                        tind,
//                        aind,
//                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX,
//                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY,
//                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta,
//                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[0].theta * 180 / PI_CONST,
//                        EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size() - 1].theta * 180 / PI_CONST,
//                        motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
//                        motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta,
//                        (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.size(),
//                        (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.size());

//             // compute linear and angular time
//             // TODO Yucheng: if necessary to calculate the cartangle?
//             // SBPL_DEBUG("Calculate linear and angular time now!");
//             double linear_distance = 0;
//             for (unsigned int i = 1; i < EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size(); i++) {
//                 // SBPL_DEBUG("point: (%f, %f, %f, %f)\n", EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].x, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].y, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].theta, EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].cartangle);
//                 double x0 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i - 1].x;
//                 double y0 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i - 1].y;
//                 double x1 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].x;
//                 double y1 = EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[i].y;
//                 double dx = x1 - x0;
//                 double dy = y1 - y0;
//                 // SBPL_DEBUG("dx, dy are %f, %f. point: (%f, %f), (%f, %f)\n", dx, dy, x0, y0, x1, y1);
//                 linear_distance += sqrt(dx * dx + dy * dy);
//                 // SBPL_DEBUG("Linear distance is %f\n", linear_distance);
//             }
//             double linear_time = linear_distance / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs;
//             // SBPL_DEBUG("Nominal linear velocity is %f\n", EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs);
//             // SBPL_DEBUG("Linear time is %f\n", linear_time);
//             double angular_distance;
//             angular_distance = fabs(computeMinUnsignedAngleDiff(
//                     DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta),
//                     DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].starttheta)));

//             double angular_time = angular_distance / ((PI_CONST / 4.0) /
//                                   EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs);
//             // SBPL_DEBUG("Angular distance is %f, Angular time is %f\n", angular_distance, angular_time);

//             // make the cost the max of the two times
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].cost =
//                     (int)(ceil(NAVXYTHETALAT_COSTMULT_MTOMM * std::max(linear_time, angular_time)));
//             // use any additional cost multiplier
//             EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;

//             // now compute the intersecting cells for this motion (including ignoring the source footprint)
//             Get2dMotionCells(
//                     EnvNAVXYTHETALATCartCfg.FootprintPolygon,
//                     EnvNAVXYTHETALATCartCfg.CartPolygon,
//                     motionprimitiveV->at(mind).intermptV,
//                     &EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV,
//                     EnvNAVXYTHETALATCartCfg.cellsize_m);

// // #if DEBUG
// //             SBPL_FPRINTF(fDeb,
// //                          "action tind=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) "
// //                          "cost=%4d (mprimID %3d: %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
// //                          tind,
// //                          aind,
// //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dX,
// //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].dY,
// //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta,
// //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[0].theta * 180 / PI_CONST,
// //                          EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV[EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intermptV.size() - 1].theta * 180 / PI_CONST, EnvNAVXYTHETALATCfg.ActionsV[tind][aind].cost,
// //                          motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
// //                          motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta,
// //                          (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].interm3DcellsV.size(),
// //                          (int)EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].intersectingcellsV.size());
// // #endif

//             // add to the list of backward actions
//             int targettheta = EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind].endtheta;
//             if (targettheta < 0) {
//                 targettheta = targettheta + EnvNAVXYTHETALATCartCfg.NumThetaDirs;
//             }
//             EnvNAVXYTHETALATCartCfg.PredActionsV[targettheta].push_back(&(EnvNAVXYTHETALATCartCfg.ActionsV[tind][aind]));
//         }


//         if (maxnumofactions < numofactions) {
//             maxnumofactions = numofactions;
//         }
    }

    // at this point we don't allow nonuniform number of actions
    if (motionprimitiveV->size() != (size_t)(EnvNAVXYTHETALATCartCfg.NumThetaDirs * CART_THETADIRS * maxnumofactions)) {
        std::stringstream ss;
        ss << "ERROR: nonuniform number of actions is not supported" <<
                " (maxnumofactions=" << maxnumofactions << " while motprims=" <<
                motionprimitiveV->size() << " thetas=" <<
                EnvNAVXYTHETALATCartCfg.NumThetaDirs;
        throw SBPL_Exception(ss.str());
    }
    
    SBPL_PRINTF("compute replanning data now!\n");

    // now compute replanning data
    ComputeReplanningDataCart();

    SBPL_PRINTF("done pre-computing action data based on motion primitives\n");
}

void EnvironmentNAVXYTHETALATTICE::Get2dMotionCells(std::vector<sbpl_2Dpt_t> polygon, std::vector<sbpl_2Dpt_t> polygon_cart, std::vector<EnvNAVXYTHETALATCART3Dpt_t> poses, std::vector<sbpl_2Dcell_t>* cells,
                         double res)
{
    // can't find any motion cells if there are no poses
    if (poses.empty()) {
        SBPL_PRINTF("Can't find any motion cells!");
        return;
    }

    //get first footprint set
    std::set<sbpl_2Dcell_t> first_cell_set;
    SBPL_PRINTF("Get first footprint set");
    CalculateFootprintForPose(polygon, polygon_cart, &first_cell_set, poses[0], res);
    // SBPL_PRINTF("The size of the first footprint set %d\n", first_cell_set.size());


    //duplicate first footprint set into motion set
    std::set<sbpl_2Dcell_t> cell_set = first_cell_set;

    //call get footprint on the rest of the points
    for (unsigned int i = 1; i < poses.size(); i++) {
        CalculateFootprintForPose(polygon, polygon_cart, &cell_set, poses[i], res);
        // SBPL_PRINTF("The size of footprint set %d\n", cell_set.size());
    }

    //convert the motion set to a vector but don't include the cells in the first footprint set
    cells->reserve(cell_set.size() - first_cell_set.size());
    for (std::set<sbpl_2Dcell_t>::iterator it = cell_set.begin(); it != cell_set.end(); it++) {
        if (first_cell_set.find(*it) == first_cell_set.end()) {
            cells->push_back(*it);
            // SBPL_PRINTF("The footprint set %f, %f\n", it->x,it->y);
        }
    }
}

bool EnvironmentNAVXYTHETALATTICE::IsValidCell(int X, int Y)
{
    return (X >= 0 && X < EnvNAVXYTHETALATCartCfg.EnvWidth_c &&
            Y >= 0 && Y < EnvNAVXYTHETALATCartCfg.EnvHeight_c &&
            EnvNAVXYTHETALATCartCfg.Grid2D[X][Y] < EnvNAVXYTHETALATCartCfg.obsthresh);
}

// trailer
bool EnvironmentNAVXYTHETALATTICE::IsWithinMapCell(int X, int Y)
{
    return (X >= 0 && X < EnvNAVXYTHETALATCartCfg.EnvWidth_c &&
            Y >= 0 && Y < EnvNAVXYTHETALATCartCfg.EnvHeight_c);
}


// trailer
bool EnvironmentNAVXYTHETALATTICE::IsValidConfiguration(int X, int Y, int Theta, int CartAngle)
{
    std::vector<sbpl_2Dcell_t> footprint;
    EnvNAVXYTHETALATCART3Dpt_t pose;

    // compute continuous pose
    pose.x = DISCXY2CONT(X, EnvNAVXYTHETALATCartCfg.cellsize_m);
    pose.y = DISCXY2CONT(Y, EnvNAVXYTHETALATCartCfg.cellsize_m);
    pose.theta = DiscTheta2ContNew(Theta);
    pose.cartangle = CartDiscTheta2Cont(CartAngle, CART_THETADIRS);

    // compute footprint cells
    CalculateFootprintForPose(
            EnvNAVXYTHETALATCartCfg.FootprintPolygon,
            EnvNAVXYTHETALATCartCfg.CartPolygon,
            &footprint,
            pose,
            EnvNAVXYTHETALATCartCfg.cellsize_m
    );

    // iterate over all footprint cells
    for (int find = 0; find < (int)footprint.size(); find++) {
        int x = footprint.at(find).x;
        int y = footprint.at(find).y;

        if (x < 0 || x >= EnvNAVXYTHETALATCartCfg.EnvWidth_c ||
            y < 0 || y >= EnvNAVXYTHETALATCartCfg.EnvHeight_c ||
            EnvNAVXYTHETALATCartCfg.Grid2D[x][y] >= EnvNAVXYTHETALATCartCfg.obsthresh)
        {
            return false;
        }
    }

    return true;
}

// trailer
void EnvironmentNAVXYTHETALATTICE::InitializeEnvConfig(std::vector<SBPL_xythetacart_mprimitive>* motionprimitiveV)
{
    // additional to configuration file initialization of EnvNAVXYTHETALATCfg if
    // necessary    

    // dXY dirs
    EnvNAVXYTHETALATCartCfg.dXY[0][0] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[0][1] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[1][0] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[1][1] = 0;
    EnvNAVXYTHETALATCartCfg.dXY[2][0] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[2][1] = 1;
    EnvNAVXYTHETALATCartCfg.dXY[3][0] = 0;
    EnvNAVXYTHETALATCartCfg.dXY[3][1] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[4][0] = 0;
    EnvNAVXYTHETALATCartCfg.dXY[4][1] = 1;
    EnvNAVXYTHETALATCartCfg.dXY[5][0] = 1;
    EnvNAVXYTHETALATCartCfg.dXY[5][1] = -1;
    EnvNAVXYTHETALATCartCfg.dXY[6][0] = 1;
    EnvNAVXYTHETALATCartCfg.dXY[6][1] = 0;
    EnvNAVXYTHETALATCartCfg.dXY[7][0] = 1;
    EnvNAVXYTHETALATCartCfg.dXY[7][1] = 1;

    EnvNAVXYTHETALATCART3Dpt_t temppose;
    temppose.x = 0.0;
    temppose.y = 0.0;
    temppose.theta = 0.0;
    temppose.cartangle = 0.0;
    std::vector<sbpl_2Dcell_t> footprint;
    // get_2d_footprint_cells(
    //         EnvNAVXYTHETALATCfg.FootprintPolygon,
    //         &footprint,
    //         temppose,
    //         EnvNAVXYTHETALATCfg.cellsize_m);
    
    // test if footprint exits
    SBPL_PRINTF("Robot footprint");
    for (size_t ii(0); ii < EnvNAVXYTHETALATCartCfg.FootprintPolygon.size(); ++ii) {
      sbpl_2Dpt_t pt;
      pt.x = EnvNAVXYTHETALATCartCfg.FootprintPolygon[ii].x;
      pt.y = EnvNAVXYTHETALATCartCfg.FootprintPolygon[ii].y;
      SBPL_PRINTF("%d: %f %f",(int)ii,pt.x,pt.y);
    }
    SBPL_PRINTF("Cart footprint");
    for (size_t ii(0); ii < EnvNAVXYTHETALATCartCfg.CartPolygon.size(); ++ii) {
      sbpl_2Dpt_t pt;
      pt.x = EnvNAVXYTHETALATCartCfg.CartPolygon[ii].x;
      pt.y = EnvNAVXYTHETALATCartCfg.CartPolygon[ii].y;
      SBPL_PRINTF("%d: %f %f",(int)ii,pt.x,pt.y);
    }

    CalculateFootprintForPose(
            EnvNAVXYTHETALATCartCfg.FootprintPolygon,
            EnvNAVXYTHETALATCartCfg.CartPolygon,
            &footprint,
            temppose,
            EnvNAVXYTHETALATCartCfg.cellsize_m);
    SBPL_PRINTF("number of cells in footprint of the robot = %d\n", (unsigned int)footprint.size());

    // for (std::vector<sbpl_2Dcell_t>::iterator it = footprint.begin(); it != footprint.end(); ++it) {
    //     SBPL_PRINTF("Footprint cell at (%d, %d)\n", it->x, it->y);
    // }

#if DEBUG
    SBPL_FPRINTF(fDeb, "footprint cells (size=%d):\n", (int)footprint.size());
    for(int i = 0; i < (int) footprint.size(); i++)
    {
        SBPL_FPRINTF(fDeb, "%d %d (cont: %.3f %.3f)\n", footprint.at(i).x, footprint.at(i).y,
                     DISCXY2CONT(footprint.at(i).x, EnvNAVXYTHETALATCartCfg.cellsize_m),
                     DISCXY2CONT(footprint.at(i).y, EnvNAVXYTHETALATCartCfg.cellsize_m));
    }
#endif

    if (motionprimitiveV == NULL) {
        SBPL_INFO("!!!!!!!!!!!!Deprecated Precompute Actions!!!!!!!!!!!!!!!!");
        // DeprecatedPrecomputeActions();
        PrecomputeActions();
    }
    else {
        PrecomputeActionswithCompleteMotionPrimitive(motionprimitiveV);
    }
}

// TODO: 31.03 not changed
// void EnvironmentNAVXYTHETACARTLATTICE::PrecomputeActions()
// {

// 	//construct list of actions
// 	ROS_DEBUG("Pre-computing action data internally using the motion primitives for a 3D kinematic planning...");
// 	EnvNAVXYTHETACARTLATCfg.ActionsV = new EnvNAVXYTHETACARTLATAction_t* [NAVXYTHETACARTLAT_THETADIRS];
// 	EnvNAVXYTHETACARTLATCfg.PredActionsV = new vector<EnvNAVXYTHETACARTLATAction_t*> [NAVXYTHETACARTLAT_THETADIRS];
// 	vector<sbpl_2Dcell_t> footprint;
// 	//iterate over source angles
// 	for(int tind = 0; tind < NAVXYTHETACARTLAT_THETADIRS; tind++)
// 	{
// 		ROS_DEBUG("processing angle %d", tind);
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind] = new EnvNAVXYTHETACARTLATAction_t[EnvNAVXYTHETACARTLATCfg.actionwidth];

// 		//compute sourcepose
// 		EnvNAVXYTHETACARTLAT3Dpt_t sourcepose;
// 		sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		sourcepose.theta = DiscTheta2Cont(tind, NAVXYTHETACARTLAT_THETADIRS);

// 		//the construction assumes that the robot first turns and then goes along this new theta
// 		int aind = 0;
// 		for(; aind < 3; aind++)
// 		{
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].starttheta = tind;
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta = (tind + aind - 1)%NAVXYTHETACARTLAT_THETADIRS; //-1,0,1
// 			double angle = DiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, NAVXYTHETACARTLAT_THETADIRS);
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX = (int)(cos(angle) + 0.5*(cos(angle)>0?1:-1));
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY = (int)(sin(angle) + 0.5*(sin(angle)>0?1:-1));
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost = (int)(ceil(NAVXYTHETACARTLAT_COSTMULT_MTOMM*EnvNAVXYTHETACARTLATCfg.cellsize_m/EnvNAVXYTHETACARTLATCfg.nominalvel_mpersecs*sqrt((double)(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX*EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX + 
// 					EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY*EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY))));

// 			//compute intersecting cells
// 			EnvNAVXYTHETACARTLAT3Dpt_t pose;
// 			pose.x = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 			pose.y = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 			pose.theta = angle;
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intermptV.clear();
// 			EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV.clear();
// 			CalculateFootprintForPose(pose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);
// 			RemoveSourceFootprint(sourcepose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);

// #if DEBUG
// 			ROS_DEBUG("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d",
// 				tind, aind, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, angle, 
// 				EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY,
// 				EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost);
// #endif

// 			//add to the list of backward actions
// 			int targettheta = EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta;
// 			if (targettheta < 0)
// 				targettheta = targettheta + NAVXYTHETACARTLAT_THETADIRS;
// 			 EnvNAVXYTHETACARTLATCfg.PredActionsV[targettheta].push_back(&(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind]));

// 		}

// 		//decrease and increase angle without movement
// 		aind = 3;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].starttheta = tind;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta = tind-1;
// 		if(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta < 0) EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta += NAVXYTHETACARTLAT_THETADIRS;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX = 0;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY = 0;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost = (int)(NAVXYTHETACARTLAT_COSTMULT_MTOMM*EnvNAVXYTHETACARTLATCfg.timetoturn45degsinplace_secs);

// 		//compute intersecting cells
// 		EnvNAVXYTHETACARTLAT3Dpt_t pose;
// 		pose.x = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		pose.y = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		pose.theta = DiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, NAVXYTHETACARTLAT_THETADIRS);
// 		pose.cartangle = CartDiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endcartangle, CART_THETADIRS);

// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intermptV.clear();
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV.clear();
// 		CalculateFootprintForPose(pose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);
// 		RemoveSourceFootprint(sourcepose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);

// #if DEBUG
// 		ROS_DEBUG("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d",
//               tind, aind, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, DiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, NAVXYTHETACARTLAT_THETADIRS),
//               EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY,
//               EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost);
// #endif

// 		//add to the list of backward actions
// 		int targettheta = EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta;
// 		if (targettheta < 0)
// 			targettheta = targettheta + NAVXYTHETACARTLAT_THETADIRS;
// 		 EnvNAVXYTHETACARTLATCfg.PredActionsV[targettheta].push_back(&(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind]));


// 		aind = 4;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].starttheta = tind;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta = (tind + 1)%NAVXYTHETACARTLAT_THETADIRS; 
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX = 0;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY = 0;
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost = (int)(NAVXYTHETACARTLAT_COSTMULT_MTOMM*EnvNAVXYTHETACARTLATCfg.timetoturn45degsinplace_secs);

// 		//compute intersecting cells
// 		pose.x = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		pose.y = DISCXY2CONT(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY, EnvNAVXYTHETACARTLATCfg.cellsize_m);
// 		pose.theta = DiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, NAVXYTHETACARTLAT_THETADIRS);
//     pose.cartangle = CartDiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endcartangle, CART_THETADIRS);
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intermptV.clear();
// 		EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV.clear();
// 		CalculateFootprintForPose(pose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);
// 		RemoveSourceFootprint(sourcepose, &EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].intersectingcellsV);


// #if DEBUG
// 		ROS_DEBUG("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d",
//               tind, aind, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, DiscTheta2Cont(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta, NAVXYTHETACARTLAT_THETADIRS),
//               EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dX, EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].dY,
//               EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].cost);
// #endif

// 		//add to the list of backward actions
// 		targettheta = EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind].endtheta;
// 		if (targettheta < 0)
// 			targettheta = targettheta + NAVXYTHETACARTLAT_THETADIRS;
// 		 EnvNAVXYTHETACARTLATCfg.PredActionsV[targettheta].push_back(&(EnvNAVXYTHETACARTLATCfg.ActionsV[tind][aind]));

// 	}

// 	//now compute replanning data
// 	ComputeReplanningDataCart();
// 	ROS_INFO("Done pre-computing action data");
// }

int EnvironmentNAVXYTHETALATTICE::GetActionCost(
    int SourceX, int SourceY, int SourceTheta, int SourceCartAngle,
    EnvNAVXYTHETALATCARTAction_t* action)
{
    sbpl_2Dcell_t cell;
    EnvNAVXYTHETALATCART3Dcell_t interm3Dcell;
    int i;

    // TODO - go over bounding box (minpt and maxpt) to test validity and skip
    // testing boundaries below, also order intersect cells so that the four
    // farthest pts go first

    if (!IsValidCell(SourceX, SourceY)) {
        return INFINITECOST;
    }
    if (!IsValidCell(SourceX + action->dX, SourceY + action->dY)) {
        return INFINITECOST;
    }

    // TODO: 1. add start and end pose collision checking for trailer?
    // 05.04 2. cost multiplicate with action length?
    if (EnvNAVXYTHETALATCartCfg.Grid2D[SourceX + action->dX][SourceY + action->dY] >=
        EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh)
    {
        return INFINITECOST;
    }

    // need to iterate over discretized center cells and compute cost based on them
    unsigned char maxcellcost = 0;
    // interm3Dcells contain the center pose of robot and the trailer
    for (i = 0; i < (int)action->interm3DcellsV.size(); i++) {
        interm3Dcell = action->interm3DcellsV.at(i);
        // SBPL_INFO("interm3DcellsV Origin %d, %d", interm3Dcell.x, interm3Dcell.y);
        interm3Dcell.x = interm3Dcell.x + SourceX;
        interm3Dcell.y = interm3Dcell.y + SourceY;
        // SBPL_INFO("interm3DcellsV in the map %d, %d", interm3Dcell.x, interm3Dcell.y);

        if (interm3Dcell.x < 0 || interm3Dcell.x >= EnvNAVXYTHETALATCartCfg.EnvWidth_c ||
            interm3Dcell.y < 0 || interm3Dcell.y >= EnvNAVXYTHETALATCartCfg.EnvHeight_c)
        {
            return INFINITECOST;
        }

        maxcellcost = __max(maxcellcost, EnvNAVXYTHETALATCartCfg.Grid2D[interm3Dcell.x][interm3Dcell.y]);

        // check that the robot is NOT in the cell at which there is no valid orientation
        if (maxcellcost >= EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh) {
            return INFINITECOST;
        }
    }

    // check collisions that for the particular footprint orientation along the action
    if (EnvNAVXYTHETALATCartCfg.FootprintPolygon.size() > 1 &&
        (int)maxcellcost >= EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh)
    {
        checks++;

        for (i = 0; i < (int)action->intersectingcellsV.size(); i++) {
            // get the cell in the map
            cell = action->intersectingcellsV.at(i);
            // SBPL_INFO("intersectingcellsV in the map %d, %d", cell.x, cell.y);

            cell.x = cell.x + SourceX;
            cell.y = cell.y + SourceY;
            // SBPL_INFO("intersectingcellsV in the map %d, %d", cell.x, cell.y);


            // check validity
            if (!IsValidCell(cell.x, cell.y)) {
                return INFINITECOST;
            }

// cost computation changed: cost = max(cost of centers of the robot along
// action) intersecting cells are only used for collision checking
//            if (EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y] > currentmaxcost) {
//              currentmaxcost = EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y];
//            }
        }
    }

    // to ensure consistency of h2D:
    maxcellcost = __max(maxcellcost, EnvNAVXYTHETALATCartCfg.Grid2D[SourceX][SourceY]);
    int currentmaxcost = (int)__max(
            maxcellcost,
            EnvNAVXYTHETALATCartCfg.Grid2D[SourceX + action->dX][SourceY + action->dY]);

    // Debug 04.04
    // SBPL_INFO("action cost with collision checking. Souce: (%d, %d, %d, %d), Cost: %d", SourceX, SourceY, SourceTheta, SourceCartAngle, action->cost * (currentmaxcost + 1));

    // use cell cost as multiplicative factor
    return action->cost * (currentmaxcost + 1);
}

double EnvironmentNAVXYTHETALATTICE::EuclideanDistance_m(int X1, int Y1, int X2, int Y2)
{
    int sqdist = ((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2));
    return EnvNAVXYTHETALATCartCfg.cellsize_m * sqrt((double)sqdist);
}

// // calculates a set of cells that correspond to the specified footprint adds
// // points to it (does not clear it beforehand)
// void EnvironmentNAVXYTHETALATTICE::CalculateFootprintForPose(
//     sbpl_xy_theta_pt_t pose, std::vector<sbpl_2Dcell_t>* footprint,
//     const std::vector<sbpl_2Dpt_t>& FootprintPolygon)
// {
//     int pind;

//     // handle special case where footprint is just a point
//     if (FootprintPolygon.size() <= 1) {
//         sbpl_2Dcell_t cell;
//         cell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETALATCfg.cellsize_m);
//         cell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETALATCfg.cellsize_m);

//         for (pind = 0; pind < (int)footprint->size(); pind++) {
//             if (cell.x == footprint->at(pind).x && cell.y == footprint->at(pind).y) break;
//         }
//         if (pind == (int)footprint->size()) {
//             footprint->push_back(cell);
//         }
//         return;
//     }

//     std::vector<sbpl_2Dpt_t> bounding_polygon;
//     unsigned int find;
//     double max_x = -INFINITECOST, min_x = INFINITECOST, max_y = -INFINITECOST, min_y = INFINITECOST;
//     sbpl_2Dpt_t pt(0, 0);
//     for (find = 0; find < FootprintPolygon.size(); find++) {
//         // rotate and translate the corner of the robot
//         pt = FootprintPolygon[find];

//         // rotate and translate the point
//         sbpl_2Dpt_t corner;
//         corner.x = cos(pose.theta) * pt.x - sin(pose.theta) * pt.y + pose.x;
//         corner.y = sin(pose.theta) * pt.x + cos(pose.theta) * pt.y + pose.y;
//         bounding_polygon.push_back(corner);
//         if (corner.x < min_x || find == 0) {
//             min_x = corner.x;
//         }
//         if (corner.x > max_x || find == 0) {
//             max_x = corner.x;
//         }
//         if (corner.y < min_y || find == 0) {
//             min_y = corner.y;
//         }
//         if (corner.y > max_y || find == 0) {
//             max_y = corner.y;
//         }
//     }

//     // initialize previous values to something that will fail the if condition
//     // during the first iteration in the for loop
//     int prev_discrete_x = CONTXY2DISC(pt.x, EnvNAVXYTHETALATCfg.cellsize_m) + 1;
//     int prev_discrete_y = CONTXY2DISC(pt.y, EnvNAVXYTHETALATCfg.cellsize_m) + 1;
//     int prev_inside = 0;
//     int discrete_x;
//     int discrete_y;

//     for (double x = min_x; x <= max_x; x += EnvNAVXYTHETALATCfg.cellsize_m / 3) {
//         for (double y = min_y; y <= max_y; y += EnvNAVXYTHETALATCfg.cellsize_m / 3) {
//             pt.x = x;
//             pt.y = y;
//             discrete_x = CONTXY2DISC(pt.x, EnvNAVXYTHETALATCfg.cellsize_m);
//             discrete_y = CONTXY2DISC(pt.y, EnvNAVXYTHETALATCfg.cellsize_m);

//             // see if we just tested this point
//             if (discrete_x != prev_discrete_x || discrete_y != prev_discrete_y || prev_inside == 0) {

//                 if (IsInsideFootprint(pt, &bounding_polygon)) {
//                     // convert to a grid point

//                     sbpl_2Dcell_t cell;
//                     cell.x = discrete_x;
//                     cell.y = discrete_y;

//                     // insert point if not there already
//                     for (pind = 0; pind < (int)footprint->size(); pind++) {
//                         if (cell.x == footprint->at(pind).x && cell.y == footprint->at(pind).y) break;
//                     }
//                     if (pind == (int)footprint->size()) {
//                         footprint->push_back(cell);
//                     }

//                     prev_inside = 1;

//                 }
//                 else {
//                     prev_inside = 0;
//                 }

//             }
//             else {

//             }

//             prev_discrete_x = discrete_x;
//             prev_discrete_y = discrete_y;
//         } // over x_min...x_max
//     }
// }

// // calculates a set of cells that correspond to the footprint of the base adds
// // points to it (does not clear it beforehand)
// void EnvironmentNAVXYTHETALATTICE::CalculateFootprintForPose(
//     sbpl_xy_theta_pt_t pose,
//     std::vector<sbpl_2Dcell_t>* footprint)
// {
//     CalculateFootprintForPose(pose, footprint, EnvNAVXYTHETALATCfg.FootprintPolygon);
// }

EnvNAVXYTHETALATCART3Dpt_t EnvironmentNAVXYTHETALATTICE::getCartCenter(EnvNAVXYTHETALATCART3Dpt_t pose, sbpl_2Dpt_t cart_center_offset)
{  
  EnvNAVXYTHETALATCART3Dpt_t cart_offset_pt;
  cart_offset_pt.x = cos(pose.theta)*EnvNAVXYTHETALATCartCfg.CartOffset.x - sin(pose.theta)*EnvNAVXYTHETALATCartCfg.CartOffset.y + pose.x;
  cart_offset_pt.y = sin(pose.theta)*EnvNAVXYTHETALATCartCfg.CartOffset.x + cos(pose.theta)*EnvNAVXYTHETALATCartCfg.CartOffset.y + pose.y;
  cart_offset_pt.theta = pose.theta + pose.cartangle;

  EnvNAVXYTHETALATCART3Dpt_t corner_cart;
  corner_cart.x = cos(cart_offset_pt.theta)*cart_center_offset.x - sin(cart_offset_pt.theta)*cart_center_offset.y + cart_offset_pt.x;
  corner_cart.y = sin(cart_offset_pt.theta)*cart_center_offset.x + cos(cart_offset_pt.theta)*cart_center_offset.y + cart_offset_pt.y;
  corner_cart.theta = cart_offset_pt.theta;
  corner_cart.cartangle = pose.cartangle;
  return corner_cart;
}

// // removes a set of cells that correspond to the specified footprint at the
// // sourcepose adds points to it (does not clear it beforehand)
// void EnvironmentNAVXYTHETALATTICE::RemoveSourceFootprint(
//     sbpl_xy_theta_pt_t sourcepose,
//     std::vector<sbpl_2Dcell_t>* footprint,
//     const std::vector<sbpl_2Dpt_t>& FootprintPolygon)
// {
//     std::vector<sbpl_2Dcell_t> sourcefootprint;

//     // compute source footprint
//     get_2d_footprint_cells(FootprintPolygon, &sourcefootprint, sourcepose, EnvNAVXYTHETALATCfg.cellsize_m);

//     // now remove the source cells from the footprint
//     for (int sind = 0; sind < (int)sourcefootprint.size(); sind++) {
//         for (int find = 0; find < (int)footprint->size(); find++) {
//             if (sourcefootprint.at(sind).x == footprint->at(find).x && sourcefootprint.at(sind).y
//                 == footprint->at(find).y) {
//                 footprint->erase(footprint->begin() + find);
//                 break;
//             }
//         } // over footprint
//     } // over source
// }

// // removes a set of cells that correspond to the footprint of the base at the
// // sourcepose adds points to it (does not clear it beforehand)
// void EnvironmentNAVXYTHETALATTICE::RemoveSourceFootprint(
//     sbpl_xy_theta_pt_t sourcepose,
//     std::vector<sbpl_2Dcell_t>* footprint)
// {
//     RemoveSourceFootprint(sourcepose, footprint, EnvNAVXYTHETALATCfg.FootprintPolygon);
// }

void EnvironmentNAVXYTHETALATTICE::EnsureHeuristicsUpdated(bool bGoalHeuristics)
{
    if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
        grid2Dsearchfromstart->search(
                EnvNAVXYTHETALATCartCfg.Grid2D,
                EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh,
                EnvNAVXYTHETALATCartCfg.StartX_c, EnvNAVXYTHETALATCartCfg.StartY_c,
                EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.EndY_c,
                SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeStartHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n", (int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.EndY_c) / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs));
    }

    if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
        grid2Dsearchfromgoal->search(
                EnvNAVXYTHETALATCartCfg.Grid2D,
                EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh,
                EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.EndY_c,
                EnvNAVXYTHETALATCartCfg.StartX_c, EnvNAVXYTHETALATCartCfg.StartY_c,
                SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeGoalHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n", (int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETALATCartCfg.StartX_c, EnvNAVXYTHETALATCartCfg.StartY_c) / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs));
    }
}

// for trailer
void EnvironmentNAVXYTHETALATTICE::ComputeHeuristicValuesCart()
{
    // whatever necessary pre-computation of heuristic values is done here
    SBPL_PRINTF("Precomputing heuristics...\n");

    // allocated 2D grid searches
    grid2Dsearchfromstart = new SBPL2DGridSearch(
            EnvNAVXYTHETALATCartCfg.EnvWidth_c, EnvNAVXYTHETALATCartCfg.EnvHeight_c,
            (float)EnvNAVXYTHETALATCartCfg.cellsize_m, blocksize, bucketsize);
    grid2Dsearchfromgoal = new SBPL2DGridSearch(
            EnvNAVXYTHETALATCartCfg.EnvWidth_c, EnvNAVXYTHETALATCartCfg.EnvHeight_c,
            (float)EnvNAVXYTHETALATCartCfg.cellsize_m, blocksize, bucketsize);

    // set OPEN type to sliding buckets
    grid2Dsearchfromstart->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
    grid2Dsearchfromgoal->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);

    SBPL_PRINTF("done\n");
}

bool EnvironmentNAVXYTHETALATTICE::CheckQuant(FILE* fOut)
{
    for (double theta = -10; theta < 10;
        theta += 2.0 * PI_CONST / EnvNAVXYTHETALATCartCfg.NumThetaDirs * 0.01)
    {
        int nTheta, nnewTheta;
        double newTheta;
        nTheta = ContTheta2DiscNew(theta);
        newTheta = DiscTheta2ContNew(nTheta);
        nnewTheta = ContTheta2DiscNew(newTheta);

        SBPL_FPRINTF(fOut, "theta=%f(%f)->%d->%f->%d\n", theta, theta * 180 / PI_CONST, nTheta, newTheta, nnewTheta);

        if (nTheta != nnewTheta) {
            SBPL_ERROR("ERROR: invalid quantization\n");
            return false;
        }
    }

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::InitializeEnv(const char* sEnvFile, const std::vector<sbpl_2Dpt_t>& perimeterptsV, const std::vector<sbpl_2Dpt_t>& cartperimeterptsV, sbpl_2Dpt_t& cart_offset, const char* sMotPrimFile)
{
	EnvNAVXYTHETALATCartCfg.FootprintPolygon = perimeterptsV;
    EnvNAVXYTHETALATCartCfg.CartPolygon = cartperimeterptsV;
    EnvNAVXYTHETALATCartCfg.CartOffset = cart_offset;

	FILE* fCfg = fopen(sEnvFile, "r");
	if(fCfg == NULL)
	{
		SBPL_ERROR("unable to open %s", sEnvFile);
        throw SBPL_Exception();
	}
	ReadConfiguration(fCfg);

	if(sMotPrimFile != NULL)
	{
		FILE* fMotPrim = fopen(sMotPrimFile, "r");
		if(fMotPrim == NULL)
		{
			SBPL_ERROR("unable to open %s", sMotPrimFile);
			throw SBPL_Exception();
		}
		if(ReadMotionPrimitivesCart(fMotPrim) == false)
		{
			SBPL_ERROR("failed to read in motion primitive file");
			throw SBPL_Exception();
		}
		InitGeneralCart(&EnvNAVXYTHETALATCartCfg.mprimV);
	}
	else
		InitGeneralCart(NULL);

	SBPL_DEBUG("size of env: %d by %d", EnvNAVXYTHETALATCartCfg.EnvWidth_c, EnvNAVXYTHETALATCartCfg.EnvHeight_c);

	return true;
}

bool EnvironmentNAVXYTHETALATTICE::InitializeEnv(const char* sEnvFile)
{
    FILE* fCfg = fopen(sEnvFile, "r");
    if (fCfg == NULL) {
        SBPL_ERROR("ERROR: unable to open %s\n", sEnvFile);
        throw SBPL_Exception();
    }
    ReadConfiguration(fCfg);
    fclose(fCfg);

    InitGeneralCart( NULL);

    return true;
}

// modified for trailer scenario
bool EnvironmentNAVXYTHETALATTICE::InitializeEnv(
    int width,
    int height,
    const unsigned char* mapdata,
    double startx, double starty, double starttheta, double startcartangle,
    double goalx, double goaly, double goaltheta, double goalcartangle,
    double goaltol_x, double goaltol_y, double goaltol_theta, double goaltol_cartangle,
    const std::vector<sbpl_2Dpt_t> & perimeterptsV,
    const sbpl_2Dpt_t & cart_offset,
    // const sbpl_2Dpt_t & cart_cp_offset,
    // const double rod_length,
    const std::vector<sbpl_2Dpt_t> & cart_perimeterptsV,
    double cellsize_m,
    double nominalvel_mpersecs,
    double timetoturn45degsinplace_secs,
    unsigned char obsthresh,
    const char* sMotPrimFile)
{
    SBPL_PRINTF("env: initialize with width=%d height=%d start=%.3f %.3f %.3f "
                "goalx=%.3f %.3f %.3f cellsize=%.3f nomvel=%.3f timetoturn=%.3f, obsthresh=%d\n",
                width, height, startx, starty, starttheta, goalx, goaly, goaltheta, cellsize_m, nominalvel_mpersecs,
                timetoturn45degsinplace_secs, obsthresh);

    SBPL_PRINTF("NOTE: goaltol parameters currently unused\n");

    SBPL_PRINTF("perimeter has size=%d\n", (unsigned int)perimeterptsV.size());

    for (int i = 0; i < (int)perimeterptsV.size(); i++) {
        SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, perimeterptsV.at(i).x, perimeterptsV.at(i).y);
    }

    SBPL_PRINTF("cart perimeter has size=%d\n", (unsigned int)cart_perimeterptsV.size());

    for (int i = 0; i < (int)cart_perimeterptsV.size(); i++) {
        SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, cart_perimeterptsV.at(i).x, cart_perimeterptsV.at(i).y);
    }

    EnvNAVXYTHETALATCartCfg.obsthresh = obsthresh;
    EnvNAVXYTHETALATCartCfg.cellsize_m = cellsize_m;
    EnvNAVXYTHETALATCartCfg.StartTheta_rad = starttheta;
    EnvNAVXYTHETALATCartCfg.EndTheta_rad = goaltheta;

    // TODO - need to set the tolerance as well

    if (sMotPrimFile != NULL) {
        FILE* fMotPrim = fopen(sMotPrimFile, "r");
        if (fMotPrim == NULL) {
            std::stringstream ss;
            ss << "ERROR: unable to open " << sMotPrimFile;
            throw SBPL_Exception(ss.str());
        }

        if (ReadMotionPrimitivesCart(fMotPrim) == false) {
            throw SBPL_Exception("ERROR: failed to read in motion primitive file");
        }
        fclose(fMotPrim);
    }
    
    // TODO: cfg for trailer
    EnvNAVXYTHETALATCartCfg.StartTheta = ContTheta2DiscNew(EnvNAVXYTHETALATCartCfg.StartTheta_rad);
    if (EnvNAVXYTHETALATCartCfg.StartTheta < 0 ||
        EnvNAVXYTHETALATCartCfg.StartTheta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs)
    {
        throw new SBPL_Exception("ERROR: illegal start coordinates for theta");
    }
    EnvNAVXYTHETALATCartCfg.EndTheta = ContTheta2DiscNew(EnvNAVXYTHETALATCartCfg.EndTheta_rad);
    if (EnvNAVXYTHETALATCartCfg.EndTheta < 0 ||
        EnvNAVXYTHETALATCartCfg.EndTheta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs)
    {
        throw new SBPL_Exception("ERROR: illegal goal coordiantes for theta");
    }


    SetConfiguration(width,
                     height,
                     mapdata,
                     CONTXY2DISC(startx, cellsize_m), 
                     CONTXY2DISC(starty, cellsize_m), 
                     EnvNAVXYTHETALATCartCfg.StartTheta,
                     ContTheta2Disc(startcartangle, CART_THETADIRS),
                     CONTXY2DISC(goalx, cellsize_m), 
                     CONTXY2DISC(goaly, cellsize_m), 
                     EnvNAVXYTHETALATCartCfg.EndTheta,
                     ContTheta2Disc(goalcartangle, CART_THETADIRS),
                     cellsize_m,
                     nominalvel_mpersecs, 
                     timetoturn45degsinplace_secs,
                     perimeterptsV,
                     cart_perimeterptsV,
                    //  cart_cp_offset,
                    //  rod_length,
                     cart_offset);

    if (EnvNAVXYTHETALATCartCfg.mprimV.size() != 0) {
        InitGeneralCart(&EnvNAVXYTHETALATCartCfg.mprimV);
    }
    else {
        InitGeneralCart( NULL);
    }

    // footprint
    PrintFootprint();

    return true;
}

//This function is inefficient and should be avoided if possible (you should
//use overloaded functions that uses a set for the cells)!
void EnvironmentNAVXYTHETALATTICE::CalculateFootprintForPose(
    std::vector<sbpl_2Dpt_t> polygon,
    std::vector<sbpl_2Dpt_t> cart_polygon,
    std::vector<sbpl_2Dcell_t>* cells, 
    EnvNAVXYTHETALATCART3Dpt_t pose,
    double res)
{
    std::set<sbpl_2Dcell_t> cell_set;
    for (unsigned int i = 0; i < cells->size(); i++)
        cell_set.insert(cells->at(i));
    CalculateFootprintForPose(polygon, cart_polygon, &cell_set, pose, res);
    cells->clear();
    cells->reserve(cell_set.size());
    //for(set<sbpl_2Dcell_t>::iterator it; it==cell_set.begin(); it++)
    //  cells->push_back(*it);
    for (std::set<sbpl_2Dcell_t>::iterator it = cell_set.begin(); it != cell_set.end(); it++) {
        cells->push_back(*it);
    }
}

// TODO Yucheng: check if the function is the same as get_2d_cells
// TODO Yucheng: change the car pushing function to tractor-trailer mode
void EnvironmentNAVXYTHETALATTICE::CalculateFootprintForPose(
    std::vector<sbpl_2Dpt_t> polygon, 
    std::vector<sbpl_2Dpt_t> cart_polygon,
    std::set<sbpl_2Dcell_t>* cells, 
    EnvNAVXYTHETALATCART3Dpt_t pose, 
    double res)
{
    // bool first_footprint;
    // if (cells->size() == 0)
    //     first_footprint = true;
    // else
    //     first_footprint = false;

    // SBPL_INFO("Size of cells at the beginning: %d, %d", cells->size(), first_footprint);

    //special case for point robot
    // SBPL_PRINTF("The size of the polygon is %d\n", (unsigned int)polygon.size());
    if (polygon.size() <= 1) {
        sbpl_2Dcell_t cell;
        cell.x = CONTXY2DISC(pose.x, res);
        cell.y = CONTXY2DISC(pose.y, res);

        cells->insert(cell);
        return;
    }

    //run bressenham line algorithm around the polygon (add them to the cells set)
    //while doing that find the min and max (x,y) and the average x and y
    double cth = cos(pose.theta);
    double sth = sin(pose.theta);

    std::vector<sbpl_2Dpt_t> disc_polygon;
    disc_polygon.reserve(polygon.size() + 1);
    int minx = INFINITECOST;
    int maxx = -INFINITECOST;
    int miny = INFINITECOST;
    int maxy = -INFINITECOST;

    //find the bounding box of the polygon
    for (unsigned int i = 0; i < polygon.size(); i++) {
        sbpl_2Dpt_t p;
        // SBPL_INFO("Polygon points: (%d, %d)", (int)(polygon[i].x > 0 ? polygon[i].x / res + 0.5 : polygon[i].x / res - 0.5), (int)(polygon[i].y > 0 ? polygon[i].y / res + 0.5 : polygon[i].y / res - 0.5));
        //		p.first = CONTXY2DISC(cth*polygon[i].x - sth*polygon[i].y + pose.x, res);
        double cx = (cth * polygon[i].x - sth * polygon[i].y + pose.x);
        double cy = (sth * polygon[i].x + cth * polygon[i].y + pose.y);
        p.x = (int)(cx > 0 ? cx / res + 0.5 : cx / res - 0.5); //(int)(cx / res + 0.5 * sign(c);
        //		p.second = CONTXY2DISC(sth*polygon[i].x + cth*polygon[i].y + pose.y, res);
        p.y = (int)(cy > 0 ? cy / res + 0.5 : cy / res - 0.5);//(int)(cy / res + 0.5);
        disc_polygon.push_back(p);
        if (p.x < minx) minx = p.x;
        if (p.x > maxx) maxx = p.x;
        if (p.y < miny) miny = p.y;
        if (p.y > maxy) maxy = p.y;
    }
    disc_polygon.push_back(disc_polygon.front());

    //make a grid big enough for the footprint
    int sizex = (maxx - minx + 1) + 2;
    int sizey = (maxy - miny + 1) + 2;
    // SBPL_PRINTF("The size of the x and y are %d / %d\n", sizex, sizey);
    int** grid = new int*[sizex];
    for (int i = 0; i < sizex; i++) {
        grid[i] = new int[sizey];
        for (int j = 0; j < sizey; j++)
            grid[i][j] = 0;
    }

    //plot line points on the grid
    for (unsigned int i = 1; i < disc_polygon.size(); i++) {
        int x0 = disc_polygon[i - 1].x - minx + 1;
        int y0 = disc_polygon[i - 1].y - miny + 1;
        int x1 = disc_polygon[i].x - minx + 1;
        int y1 = disc_polygon[i].y - miny + 1;

        //bressenham (add the line cells to the set and to a vector)
        bool steep = abs(y1 - y0) > abs(x1 - x0);
        if (steep) {
            int temp = x0;
            x0 = y0;
            y0 = temp;
            temp = x1;
            x1 = y1;
            y1 = temp;
        }
        if (x0 > x1) {
            int temp = x0;
            x0 = x1;
            x1 = temp;
            temp = y0;
            y0 = y1;
            y1 = temp;
        }
        int deltax = x1 - x0;
        int deltay = abs(y1 - y0);
        int error = deltax / 2;
        int ystep = (y0 < y1 ? 1 : -1);
        int y = y0;
        for (int x = x0; x <= x1; x++) {
            if (steep) {
                grid[y][x] = 1;
                cells->insert(sbpl_2Dcell_t(y - 1 + minx, x - 1 + miny));
            }
            else {
                grid[x][y] = 1;
                cells->insert(sbpl_2Dcell_t(x - 1 + minx, y - 1 + miny));
            }
            int last_error = error;
            error -= deltay;
            if (error < 0 && x != x1) {
                //make sure we can't have a diagonal line (the 8-connected bfs will leak through)

                int tempy = y;
                int tempx = x;
                if (last_error < -error)
                    tempy += ystep;
                else
                    tempx += 1;
                if (steep) {
                    grid[tempy][tempx] = 1;
                    cells->insert(sbpl_2Dcell_t(tempy - 1 + minx, tempx - 1 + miny));
                }
                else {
                    grid[tempx][tempy] = 1;
                    cells->insert(sbpl_2Dcell_t(tempx - 1 + minx, tempy - 1 + miny));
                }

                y += ystep;
                error += deltax;
            }
        }
    }

    //run a 2d bfs from the average (x,y)
    sbpl_bfs_2d bfs(sizex, sizey, 1);
    bfs.compute_distance_from_point(grid, 0, 0);

    for (int i = 0; i < sizex; i++)
        delete[] grid[i];
    delete[] grid;

    //add all cells expanded to the cells set
    for (int i = 1; i < sizex - 1; i++) {
        for (int j = 1; j < sizey - 1; j++) {
            if (bfs.get_distance(i, j) < 0) cells->insert(sbpl_2Dcell_t(i - 1 + minx, j - 1 + miny));
        }
    }
    ///////////////////////////
    // Do the same for the cart
    ///////////////////////////
    std::vector<sbpl_2Dpt_t> disc_polygon_cart;
    disc_polygon_cart.reserve(cart_polygon.size() + 1);
    int minx_cart = INFINITECOST;
    int maxx_cart = -INFINITECOST;
    int miny_cart = INFINITECOST;
    int maxy_cart = -INFINITECOST;

    sbpl_2Dpt_t cart_offset_pt;
    cart_offset_pt.x = cth*EnvNAVXYTHETALATCartCfg.CartOffset.x - sth*EnvNAVXYTHETALATCartCfg.CartOffset.y + pose.x;
    cart_offset_pt.y = sth*EnvNAVXYTHETALATCartCfg.CartOffset.x + cth*EnvNAVXYTHETALATCartCfg.CartOffset.y + pose.y;
    double cart_angle_global = pose.theta + pose.cartangle;

    // double cth_cart = cos(pose.theta);
    // double sth_cart = sin(pose.theta);
    double cth_cart = cos(cart_angle_global);
    double sth_cart = sin(cart_angle_global);

    //find the bounding box of the polygon
    for (unsigned int i = 0; i < cart_polygon.size(); i++) {
        sbpl_2Dpt_t p;
        // SBPL_INFO("Cart Polygon points: (%d, %d)", (int)(cart_polygon[i].x > 0 ? cart_polygon[i].x / res + 0.5 : cart_polygon[i].x / res - 0.5), (int)(cart_polygon[i].y > 0 ? cart_polygon[i].y / res + 0.5 : cart_polygon[i].y / res - 0.5));

        //		p.first = CONTXY2DISC(cth*polygon[i].x - sth*polygon[i].y + pose.x, res);
        double cx = (cth_cart * cart_polygon[i].x - sth_cart * cart_polygon[i].y + cart_offset_pt.x);
        double cy = (sth_cart * cart_polygon[i].x + cth_cart * cart_polygon[i].y + cart_offset_pt.y);
        p.x = (int)(cx > 0 ? cx / res + 0.5 : cx / res - 0.5); //(int)(cx / res + 0.5 * sign(c);
        //		p.second = CONTXY2DISC(sth*polygon[i].x + cth*polygon[i].y + pose.y, res);
        p.y = (int)(cy > 0 ? cy / res + 0.5 : cy / res - 0.5);//(int)(cy / res + 0.5);
        disc_polygon_cart.push_back(p);
        if (p.x < minx_cart) minx_cart = p.x;
        if (p.x > maxx_cart) maxx_cart = p.x;
        if (p.y < miny_cart) miny_cart = p.y;
        if (p.y > maxy_cart) maxy_cart = p.y;
    }
    disc_polygon_cart.push_back(disc_polygon_cart.front());

    //make a grid big enough for the footprint
    int sizex_cart = (maxx_cart - minx_cart + 1) + 2;
    int sizey_cart = (maxy_cart - miny_cart + 1) + 2;
    // SBPL_PRINTF("The cart size of the x and y are %d / %d\n", sizex_cart, sizey_cart);
    int** grid_cart = new int*[sizex_cart];
    for (int i = 0; i < sizex_cart; i++) {
        grid_cart[i] = new int[sizey_cart];
        for (int j = 0; j < sizey_cart; j++)
            grid_cart[i][j] = 0;
    }

    //plot line points on the grid
    for (unsigned int i = 1; i < disc_polygon_cart.size(); i++) {
        int x0 = disc_polygon_cart[i - 1].x - minx_cart + 1;
        int y0 = disc_polygon_cart[i - 1].y - miny_cart + 1;
        int x1 = disc_polygon_cart[i].x - minx_cart + 1;
        int y1 = disc_polygon_cart[i].y - miny_cart + 1;

        //bressenham (add the line cells to the set and to a vector)
        bool steep = abs(y1 - y0) > abs(x1 - x0);
        if (steep) {
            int temp = x0;
            x0 = y0;
            y0 = temp;
            temp = x1;
            x1 = y1;
            y1 = temp;
        }
        if (x0 > x1) {
            int temp = x0;
            x0 = x1;
            x1 = temp;
            temp = y0;
            y0 = y1;
            y1 = temp;
        }
        int deltax = x1 - x0;
        int deltay = abs(y1 - y0);
        int error = deltax / 2;
        int ystep = (y0 < y1 ? 1 : -1);
        int y = y0;
        for (int x = x0; x <= x1; x++) {
            if (steep) {
                grid_cart[y][x] = 1;
                cells->insert(sbpl_2Dcell_t(y - 1 + minx_cart, x - 1 + miny_cart));
            }
            else {
                grid_cart[x][y] = 1;
                cells->insert(sbpl_2Dcell_t(x - 1 + minx_cart, y - 1 + miny_cart));
            }
            int last_error = error;
            error -= deltay;
            if (error < 0 && x != x1) {
                //make sure we can't have a diagonal line (the 8-connected bfs will leak through)

                int tempy = y;
                int tempx = x;
                if (last_error < -error)
                    tempy += ystep;
                else
                    tempx += 1;
                if (steep) {
                    grid_cart[tempy][tempx] = 1;
                    cells->insert(sbpl_2Dcell_t(tempy - 1 + minx_cart, tempx - 1 + miny_cart));
                }
                else {
                    grid_cart[tempx][tempy] = 1;
                    cells->insert(sbpl_2Dcell_t(tempx - 1 + minx_cart, tempy - 1 + miny_cart));
                }

                y += ystep;
                error += deltax;
            }
        }
    }

    //run a 2d bfs from the average (x,y)
    sbpl_bfs_2d bfs_cart(sizex_cart, sizey_cart, 1);
    bfs_cart.compute_distance_from_point(grid_cart, 0, 0);

    for (int i = 0; i < sizex_cart; i++)
        delete[] grid_cart[i];
    delete[] grid_cart;

    //add all cells expanded to the cells set
    for (int i = 1; i < sizex_cart - 1; i++) {
        for (int j = 1; j < sizey_cart - 1; j++) {
            if (bfs_cart.get_distance(i, j) < 0) cells->insert(sbpl_2Dcell_t(i - 1 + minx_cart, j - 1 + miny_cart));
        }
    }

    // Debug: if all footprint is calculated correctly
    // std::ofstream outfile;
    // outfile.open ("/home/robot/slamdog/pipeline_ws/src/nav_potential_ws/src/data.csv", std::ofstream::trunc);
    // outfile << "The pose of the robot: ,"<<pose.x <<", "<< pose.y << ", "<< pose.theta << ", "<<  pose.cartangle << ", " << minx << ","<< maxx <<", "<< miny<<", "<< maxy <<std::endl;
    // for ( const sbpl_2Dcell_t &p : *cells )
    // {
    //     // SBPL_INFO("the footprint cell: (%d, %d)", p.x, p.y);
    //     outfile << "The footprint cell: ,"<< p.x <<", "<< p.y << std::endl;
    // }
    // outfile.close();
    // if(pose.theta>1.5 && pose.theta<3 && first_footprint)
    //     throw SBPL_Exception("ERROR: stop for test");
    
    // first_footprint = false;


    // outFile.open("data.csv",std::ios::out);
    // outFile.open("/home/robot/slamdog/pipeline_ws/src/nav_potential_ws/src/data.csv", std::ofstream::out);
    // SBPL_INFO("The pose of robot: (%f, %f, %f, %f), the size of footprint: %d\n",pose.x, pose.y, pose.theta, pose.cartangle, cells->size());
    // // outFile << "The pose of the robot: "<<pose.x <<", "<< pose.y << ", "<< pose.theta << ", "<<  pose.cartangle << std::endl;
    // for ( const sbpl_2Dcell_t &p : *cells )
    // {
    //     SBPL_INFO("the footprint cell: (%d, %d)", p.x, p.y);
    //     //outFile << "The footprint cell: ("<< p.x <<", "<< p.y << std::endl;
    // }
    // outFile.close();
}


void EnvironmentNAVXYTHETALATTICE::PrintFootprint()
{
  double cart_angle;
  for(unsigned int i=0; i < CART_THETADIRS; i++)
  {
    cart_angle = CartDiscTheta2Cont(i, CART_THETADIRS);
    SBPL_PRINTF("Cart angle: discretized: %d, actual: %f",i,cart_angle);
  }

  int cart_angle_d;
  for(int i=-10; i <= 10; i++)
  {
    cart_angle_d = CartContTheta2Disc(i*MAX_CART_ANGLE/10.0, CART_THETADIRS);
    SBPL_PRINTF("Cart angle: actual: %f, discretized: %d",i*MAX_CART_ANGLE/10.0,cart_angle_d);
  }
}

// trailer
bool EnvironmentNAVXYTHETALATTICE::InitGeneralCart(
    std::vector<SBPL_xythetacart_mprimitive>* motionprimitiveV)
{
    // Initialize other parameters of the environment
    InitializeEnvConfig(motionprimitiveV);

    // initialize Environment
    InitializeEnvironmentCart();

    // pre-compute heuristics
    ComputeHeuristicValuesCart();

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::InitializeMDPCfg(MDPConfig* MDPCfg)
{
    // initialize MDPCfg with the start and goal ids
    MDPCfg->goalstateid = EnvNAVXYTHETALATCART.goalstateid;
    MDPCfg->startstateid = EnvNAVXYTHETALATCART.startstateid;

    return true;
}

void EnvironmentNAVXYTHETALATTICE::PrintHeuristicValues()
{
    const char* heur = "heur.txt";
    FILE* fHeur = SBPL_FOPEN(heur, "w");
    if (fHeur == NULL) {
        throw SBPL_Exception("ERROR: could not open debug file to write heuristic");
    }
    SBPL2DGridSearch* grid2Dsearch = NULL;

    for (int i = 0; i < 2; i++) {
        if (i == 0 && grid2Dsearchfromstart != NULL) {
            grid2Dsearch = grid2Dsearchfromstart;
            SBPL_FPRINTF(fHeur, "start heuristics:\n");
        }
        else if (i == 1 && grid2Dsearchfromgoal != NULL) {
            grid2Dsearch = grid2Dsearchfromgoal;
            SBPL_FPRINTF(fHeur, "goal heuristics:\n");
        }
        else {
            continue;
        }

        for (int y = 0; y < EnvNAVXYTHETALATCartCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvNAVXYTHETALATCartCfg.EnvWidth_c; x++) {
                if (grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y) < INFINITECOST) {
                    SBPL_FPRINTF(fHeur, "%5d ", grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y));
                }
                else {
                    SBPL_FPRINTF(fHeur, "XXXXX ");
                }
            }
            SBPL_FPRINTF(fHeur, "\n");
        }
    }
    SBPL_FCLOSE(fHeur);
}

void EnvironmentNAVXYTHETALATTICE::SetAllPreds(CMDPSTATE* state)
{
    // implement this if the planner needs access to predecessors
    throw SBPL_Exception("ERROR in EnvNAVXYTHETALAT... function: SetAllPreds is undefined");
}

void EnvironmentNAVXYTHETALATTICE::GetSuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV)
{
    GetSuccs(SourceStateID, SuccIDV, CostV, NULL);
}
void EnvironmentNAVXYTHETALATTICE::GetLazySuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazySuccs(SourceStateID, SuccIDV, CostV, isTrueCost, NULL);
}

void EnvironmentNAVXYTHETALATTICE::GetSuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV)
{
    GetSuccsWithUniqueIds(SourceStateID, SuccIDV, CostV, NULL);
}

void EnvironmentNAVXYTHETALATTICE::GetLazySuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazySuccsWithUniqueIds(SourceStateID, SuccIDV, CostV, isTrueCost, NULL);
}

const EnvNAVXYTHETALATConfig_Cart_t* EnvironmentNAVXYTHETALATTICE::GetEnvNavConfig()
{
    return &EnvNAVXYTHETALATCartCfg;
}

bool EnvironmentNAVXYTHETALATTICE::UpdateCost(
    int x,
    int y,
    unsigned char newcost)
{
    EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = newcost;

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::UpdateCostCart(
    int x,
    int y,
    unsigned char newcost)
{
    // SBPL_INFO("Cost updated for cell %d %d from old cost=%d to new cost=%u", x,y,EnvNAVXYTHETALATCartCfg.Grid2D[x][y], newcost);

    EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = newcost;

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::UpdateCostCart(
    unsigned char newcost)
{
    SBPL_INFO("Cost updated for cell from old cost to new cost=%u", newcost);

    // EnvNAVXYTHETALATCartCfg.Grid2D[x][y] = newcost;

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

bool EnvironmentNAVXYTHETALATTICE::SetMap(const unsigned char* mapdata)
{
    int xind = -1, yind = -1;

    for (xind = 0; xind < EnvNAVXYTHETALATCartCfg.EnvWidth_c; xind++) {
        for (yind = 0; yind < EnvNAVXYTHETALATCartCfg.EnvHeight_c; yind++) {
            EnvNAVXYTHETALATCartCfg.Grid2D[xind][yind] = mapdata[xind + yind * EnvNAVXYTHETALATCartCfg.EnvWidth_c];
        }
    }

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

void EnvironmentNAVXYTHETALATTICE::PrintEnv_Config(FILE* fOut)
{
    // implement this if the planner needs to print out EnvNAVXYTHETALAT. configuration
    throw SBPL_Exception("ERROR in EnvNAVXYTHETALAT... function: PrintEnv_Config is undefined");
}

void EnvironmentNAVXYTHETALATTICE::Set2DBlockSize(int BlockSize)
{
    blocksize = BlockSize;
}

void EnvironmentNAVXYTHETALATTICE::Set2DBucketSize(int BucketSize)
{
    bucketsize = BucketSize;
}

void EnvironmentNAVXYTHETALATTICE::PrintTimeStat(FILE* fOut)
{
#if TIME_DEBUG
    SBPL_FPRINTF(fOut, "time3_addallout = %f secs, time_gethash = %f secs, time_createhash = %f secs, "
                 "time_getsuccs = %f\n",
                 time3_addallout/(double)CLOCKS_PER_SEC, time_gethash/(double)CLOCKS_PER_SEC,
                 time_createhash/(double)CLOCKS_PER_SEC, time_getsuccs/(double)CLOCKS_PER_SEC);
#endif
}

bool EnvironmentNAVXYTHETALATTICE::IsObstacle(int x, int y)
{
#if DEBUG
    SBPL_FPRINTF(fDeb, "Status of cell %d %d is queried. Its cost=%d\n", x,y,EnvNAVXYTHETALATCartCfg.Grid2D[x][y]);
#endif

    return EnvNAVXYTHETALATCartCfg.Grid2D[x][y] >= EnvNAVXYTHETALATCartCfg.obsthresh;
}

void EnvironmentNAVXYTHETALATTICE::GetEnvParms(
    int *size_x, int *size_y, int* num_thetas,
    double* startx, double* starty, double* starttheta,
    double* goalx, double* goaly, double* goaltheta,
    double* cellsize_m,
    double* nominalvel_mpersecs,
    double* timetoturn45degsinplace_secs,
    unsigned char* obsthresh,
    std::vector<SBPL_xythetacart_mprimitive>* mprimitiveV)
{
    *num_thetas = EnvNAVXYTHETALATCartCfg.NumThetaDirs;
    GetEnvParms(
            size_x, size_y,
            startx, starty, starttheta,
            goalx, goaly, goaltheta,
            cellsize_m,
            nominalvel_mpersecs, timetoturn45degsinplace_secs,
            obsthresh,
            mprimitiveV);
}

void EnvironmentNAVXYTHETALATTICE::GetEnvParms(
    int* size_x, int* size_y,
    double* startx, double* starty, double* starttheta,
    double* goalx, double* goaly, double* goaltheta,
    double* cellsize_m,
    double* nominalvel_mpersecs, double* timetoturn45degsinplace_secs,
    unsigned char* obsthresh,
    std::vector<SBPL_xythetacart_mprimitive>* mprimitiveV)
{
    *size_x = EnvNAVXYTHETALATCartCfg.EnvWidth_c;
    *size_y = EnvNAVXYTHETALATCartCfg.EnvHeight_c;

    *startx = DISCXY2CONT(EnvNAVXYTHETALATCartCfg.StartX_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
    *starty = DISCXY2CONT(EnvNAVXYTHETALATCartCfg.StartY_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
    *starttheta = DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.StartTheta);
    *goalx = DISCXY2CONT(EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
    *goaly = DISCXY2CONT(EnvNAVXYTHETALATCartCfg.EndY_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
    *goaltheta = DiscTheta2ContNew(EnvNAVXYTHETALATCartCfg.EndTheta);

    *cellsize_m = EnvNAVXYTHETALATCartCfg.cellsize_m;
    *nominalvel_mpersecs = EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs;
    *timetoturn45degsinplace_secs = EnvNAVXYTHETALATCartCfg.timetoturn45degsinplace_secs;

    *obsthresh = EnvNAVXYTHETALATCartCfg.obsthresh;

    *mprimitiveV = EnvNAVXYTHETALATCartCfg.mprimV;
}

bool EnvironmentNAVXYTHETALATTICE::PoseContToDisc(
    double px, double py, double pth, int &ix, int &iy, int &ith) const
{
    ix = CONTXY2DISC(px, EnvNAVXYTHETALATCartCfg.cellsize_m);
    iy = CONTXY2DISC(py, EnvNAVXYTHETALATCartCfg.cellsize_m);
    ith = ContTheta2DiscNew(pth);
    return (pth >= -2 * PI_CONST) && (pth <= 2 * PI_CONST) &&
            (ix >= 0) && (ix < EnvNAVXYTHETALATCartCfg.EnvWidth_c) &&
            (iy >= 0) && (iy < EnvNAVXYTHETALATCartCfg.EnvHeight_c);
}

bool EnvironmentNAVXYTHETALATTICE::PoseDiscToCont(
    int ix, int iy, int ith, double &px, double &py, double &pth) const
{
    px = DISCXY2CONT(ix, EnvNAVXYTHETALATCartCfg.cellsize_m);
    py = DISCXY2CONT(iy, EnvNAVXYTHETALATCartCfg.cellsize_m);
    pth = normalizeAngle(DiscTheta2ContNew(ith));
    return (ith >= 0) && (ith < EnvNAVXYTHETALATCartCfg.NumThetaDirs) &&
            (ix >= 0) && (ix < EnvNAVXYTHETALATCartCfg.EnvWidth_c) &&
            (iy >= 0) && (iy < EnvNAVXYTHETALATCartCfg.EnvHeight_c);
}

// unsigned char EnvironmentNAVXYTHETALATTICE::GetMapCost(int x, int y)
// {
//     return EnvNAVXYTHETALATCfg.Grid2D[x][y];
// }

unsigned char EnvironmentNAVXYTHETALATTICE::GetMapCost(int x, int y)
{
    return EnvNAVXYTHETALATCartCfg.Grid2D[x][y];
}

bool EnvironmentNAVXYTHETALATTICE::SetEnvParameter(
    const char* parameter,
    int value)
{
    if (EnvNAVXYTHETALATCART.bInitialized) {
        SBPL_ERROR("ERROR: all parameters must be set before initialization of the environment\n");
        return false;
    }

    SBPL_PRINTF("setting parameter %s to %d\n", parameter, value);

    if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh = (unsigned char)value;
    }
    else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh = value;
    }
    else if (strcmp(parameter, "cost_obsthresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvNAVXYTHETALATCartCfg.obsthresh = (unsigned char)value;
    }
    else {
        SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
        return false;
    }

    return true;
}

int EnvironmentNAVXYTHETALATTICE::GetEnvParameter(const char* parameter)
{
    if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
        return (int)EnvNAVXYTHETALATCartCfg.cost_inscribed_thresh;
    }
    else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
        return (int)EnvNAVXYTHETALATCartCfg.cost_possibly_circumscribed_thresh;
    }
    else if (strcmp(parameter, "cost_obsthresh") == 0) {
        return (int)EnvNAVXYTHETALATCartCfg.obsthresh;
    }
    else {
        std::stringstream ss;
        ss << "ERROR: invalid parameter " << parameter;
        throw SBPL_Exception(ss.str());
    }
}

EnvironmentNAVXYTHETALAT::~EnvironmentNAVXYTHETALAT()
{
    SBPL_PRINTF("destroying XYTHETALAT\n");

    // delete the states themselves first
    // for (int i = 0; i < (int)StateID2CoordTable.size(); i++) {
    //     delete StateID2CoordTable.at(i);
    //     StateID2CoordTable.at(i) = NULL;
    // }
    // StateID2CoordTable.clear();

    // // delete hashtable
    // if (Coord2StateIDHashTable != NULL) {
    //     delete[] Coord2StateIDHashTable;
    //     Coord2StateIDHashTable = NULL;
    // }
    // if (Coord2StateIDHashTable_lookup != NULL) {
    //     delete[] Coord2StateIDHashTable_lookup;
    //     Coord2StateIDHashTable_lookup = NULL;
    // }

    // for trailor
    for (int i = 0; i < (int)StateID2CoordTableCart.size(); i++) {
        delete StateID2CoordTableCart.at(i);
        StateID2CoordTableCart.at(i) = NULL;
    }
    StateID2CoordTableCart.clear();

    // delete hashtable
    if (Coord2StateIDHashTableCart != NULL) {
        delete[] Coord2StateIDHashTableCart;
        Coord2StateIDHashTableCart = NULL;
    }
    if (Coord2StateIDHashTableCart_lookup != NULL) {
        delete[] Coord2StateIDHashTableCart_lookup;
        Coord2StateIDHashTableCart_lookup = NULL;
    }
}

void EnvironmentNAVXYTHETALAT::GetCoordFromState(
    int stateID, int& x, int& y, int& theta, int &cartangle) const
{
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[stateID];
    x = HashEntry->X;
    y = HashEntry->Y;
    theta = HashEntry->Theta;
    cartangle = HashEntry->CartAngle;
    SBPL_PRINTF("The stateID %d corresponding pose is %d, %d, %d, %d.\n", stateID, HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle);
}

// void EnvironmentNAVXYTHETALAT::GetCoordFromState(
//     int stateID, int& x, int& y, int& theta) const
// {
//     EnvNAVXYTHETALATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
//     x = HashEntry->X;
//     y = HashEntry->Y;
//     theta = HashEntry->Theta;
// }

int EnvironmentNAVXYTHETALAT::GetStateFromCoord(int x, int y, int theta, int cartangle)
{
    EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
    if ((OutHashEntry = (this->*GetHashEntryCart)(x, y, theta, cartangle)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntryCart)(x, y, theta, cartangle);
    }
    return OutHashEntry->stateID;
}

// not changed for trailer config yet
void EnvironmentNAVXYTHETALAT::GetActionsFromStateIDPath(
    std::vector<int>* stateIDPath,
    std::vector<EnvNAVXYTHETALATAction_t>* action_list)
{
    std::vector<EnvNAVXYTHETALATAction_t*> actionV;
    std::vector<int> CostV;
    std::vector<int> SuccIDV;
    int targetx_c, targety_c, targettheta_c, targetphi_c;
    int sourcex_c, sourcey_c, sourcetheta_c, sourcephi_c;

    SBPL_PRINTF("checks=%ld\n", checks);

    action_list->clear();

    for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
        int sourceID = stateIDPath->at(pind);
        int targetID = stateIDPath->at(pind + 1);

        // get successors and pick the target via the cheapest action
        SuccIDV.clear();
        CostV.clear();
        actionV.clear();
        // GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

        int bestcost = INFINITECOST;
        int bestsind = -1;

        for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
            if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
                bestcost = CostV[sind];
                bestsind = sind;
            }
        }
        if (bestsind == -1) {
            SBPL_ERROR("ERROR: successor not found for transition");
            GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcephi_c);
            GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, sourcephi_c);
            SBPL_PRINTF("%d %d %d -> %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, sourcephi_c, targetx_c, targety_c, targettheta_c, targetphi_c);
            throw SBPL_Exception("ERROR: successor not found for transition");
        }

#if DEBUG
        SBPL_FPRINTF(fDeb, "Start: %.3f %.3f %.3f Target: %.3f %.3f %.3f Prim ID, Start Theta: %d %d\n",
            sourcex_c, sourcey_c, sourcetheta_c, sourcephi_c,
            targetx_c, targety_c, targettheta_c, targetphi_c,
            actionV[bestsind]->aind, actionV[bestsind]->starttheta);
#endif

        action_list->push_back(*(actionV[bestsind]));
    }
}

// for trailer
void EnvironmentNAVXYTHETALAT::ConvertStateIDPathintoXYThetaPath(
    std::vector<int>* stateIDPath,
    std::vector<EnvNAVXYTHETALATCART3Dpt_t>* xythetaPath)
{
    std::vector<EnvNAVXYTHETALATCARTAction_t*> actionV;
    std::vector<int> CostV;
    std::vector<int> SuccIDV;
    int targetx_c, targety_c, targettheta_c, targetcartangle_c;
    int sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c;

    SBPL_PRINTF("checks=%ld\n", checks);

    xythetaPath->clear();

// #if DEBUG
//     SBPL_FPRINTF(fDeb, "converting stateid path into coordinates:\n");
// #endif

    for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
        int sourceID = stateIDPath->at(pind);
        int targetID = stateIDPath->at(pind + 1);

// #if DEBUG
//         GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c);
// #endif

        // get successors and pick the target via the cheapest action
        SuccIDV.clear();
        CostV.clear();
        actionV.clear();
        GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

        int bestcost = INFINITECOST;
        int bestsind = -1;

// #if DEBUG
        GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c);
        GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetcartangle_c);
        SBPL_INFO("looking for %d %d %d %d -> %d %d %d %d(numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c, targetx_c, targety_c, targettheta_c, targetcartangle_c, (int)SuccIDV.size());
        // SBPL_FPRINTF(fDeb, "looking for %d %d %d %d -> %d %d %d %d(numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c, targetx_c, targety_c, targettheta_c, targetcartangle_c, (int)SuccIDV.size());
// #endif

        for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
// #if DEBUG
            // int x_c, y_c, theta_c, cartangle_c;
            // GetCoordFromState(SuccIDV[sind], x_c, y_c, theta_c, cartangle_c);
            // SBPL_FPRINTF(fDeb, "succ: %d %d %d %d\n", x_c, y_c, theta_c, cartangle_c);
// #endif
            if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
                bestcost = CostV[sind];
                bestsind = sind;
            }
        }
        if (bestsind == -1) {
            SBPL_ERROR("ERROR: successor not found for transition");
            GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c);
            GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetcartangle_c);
            SBPL_PRINTF("%d %d %d %d-> %d %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c, targetx_c, targety_c, targettheta_c, targetcartangle_c);
            throw SBPL_Exception("ERROR: successor not found for transition");
        }

        // now push in the actual path
        GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcecartangle_c);
        double sourcex, sourcey;
        sourcex = DISCXY2CONT(sourcex_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
        sourcey = DISCXY2CONT(sourcey_c, EnvNAVXYTHETALATCartCfg.cellsize_m);
        // TODO - when there are no motion primitives we should still print source state
        for (int ipind = 0; ipind < ((int)actionV[bestsind]->intermptV.size()) - 1; ipind++) {
            // translate appropriately
            EnvNAVXYTHETALATCART3Dpt_t intermpt = actionV[bestsind]->intermptV[ipind];
            intermpt.x += sourcex;
            intermpt.y += sourcey;
        // SBPL_PRINTF("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!intermpt is %f, %f\n", intermpt.x, intermpt.x);
        // only for debug
        GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetcartangle_c);


// #if DEBUG
//             int nx = CONTXY2DISC(intermpt.x, EnvNAVXYTHETALATCfg.cellsize_m);
//             int ny = CONTXY2DISC(intermpt.y, EnvNAVXYTHETALATCfg.cellsize_m);
//             int ntheta;
//             ntheta = ContTheta2DiscNew(intermpt.theta);

//             SBPL_FPRINTF(fDeb, "%.3f %.3f %.3f (%d %d %d cost=%d) ", intermpt.x, intermpt.y, intermpt.theta, nx, ny, ntheta,  EnvNAVXYTHETALATCfg.Grid2D[nx][ny]);

//             if (ipind == 0) {
//                 SBPL_FPRINTF(fDeb, "first (heur=%d)\n", GetStartHeuristic(sourceID));
//             }
//             else {
//                 SBPL_FPRINTF(fDeb, "\n");
//             }
// #endif
            // store
            xythetaPath->push_back(intermpt);
        }
    }
}

// trailer
int EnvironmentNAVXYTHETALAT::SetGoal(double x_m, double y_m, double theta_rad, double cartangle_rad)
{
    int x = CONTXY2DISC(x_m, EnvNAVXYTHETALATCartCfg.cellsize_m);
    int y = CONTXY2DISC(y_m, EnvNAVXYTHETALATCartCfg.cellsize_m);
    int theta = ContTheta2DiscNew(theta_rad);
    int cartangle = CartContTheta2Disc(cartangle_rad, CART_THETADIRS);

    SBPL_PRINTF("env: setting goal to %.3f %.3f %.3f %.3f (%d %d %d %d)\n", x_m, y_m, theta_rad, cartangle_rad, x, y, theta, cartangle);

    if (!IsWithinMapCell(x, y)) {
        SBPL_ERROR("ERROR: trying to set a goal cell %d %d that is outside of map\n", x, y);
        return -1;
    }

    if (!IsValidConfiguration(x, y, theta, cartangle)) {
        SBPL_PRINTF("WARNING: goal configuration is invalid\n");
    }

    EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
    if ((OutHashEntry = (this->*GetHashEntryCart)(x, y, theta, cartangle)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntryCart)(x, y, theta, cartangle);
    }

    // need to recompute start heuristics?
    if (EnvNAVXYTHETALATCART.goalstateid != OutHashEntry->stateID) {
        // because termination condition may not plan all the way to the new goal
        bNeedtoRecomputeStartHeuristics = true;

        // because goal heuristics change
        bNeedtoRecomputeGoalHeuristics = true;
    }

    EnvNAVXYTHETALATCART.goalstateid = OutHashEntry->stateID;

    EnvNAVXYTHETALATCartCfg.EndX_c = x;
    EnvNAVXYTHETALATCartCfg.EndY_c = y;
    EnvNAVXYTHETALATCartCfg.EndTheta = theta;
    EnvNAVXYTHETALATCartCfg.EndCartAngle = cartangle;

    // SBPL_INFO("env: setting goal ID %d (%d %d %d %d)\n", EnvNAVXYTHETALATCART.goalstateid, EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.EndY_c, EnvNAVXYTHETALATCartCfg.EndTheta, EnvNAVXYTHETALATCartCfg.EndCartAngle);

    return EnvNAVXYTHETALATCART.goalstateid;
}

int EnvironmentNAVXYTHETALAT::SetStart(double x_m, double y_m, double theta_rad, double cartangle_rad)
{
    // SBPL_INFO("setting start! cellsize_m %d", EnvNAVXYTHETALATCartCfg.cellsize_m);
    int x = CONTXY2DISC(x_m, EnvNAVXYTHETALATCartCfg.cellsize_m);
    int y = CONTXY2DISC(y_m, EnvNAVXYTHETALATCartCfg.cellsize_m);

    int theta = ContTheta2DiscNew(theta_rad);
    int cartangle = CartContTheta2Disc(cartangle_rad, CART_THETADIRS);

    if (!IsWithinMapCell(x, y)) {
        SBPL_ERROR("ERROR: trying to set a start cell %d %d that is outside of map\n", x, y);
        return -1;
    }

    SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f (%d %d %d %d)\n", x_m, y_m, theta_rad, cartangle_rad, x, y, theta, cartangle);

    // SBPL_PRINTF("is valid configuration: %d\n", IsValidConfiguration(x, y, theta, cartangle));

    if (!IsValidConfiguration(x, y, theta, cartangle)) {
        SBPL_PRINTF("WARNING: start configuration %d %d %d %d is invalid\n", x, y, theta, cartangle);
    }

    EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
    // TODO Yucheng: find why the function cannot be called!!!
    // OutHashEntry = (this->*GetHashEntryCart)(x, y, theta, cartangle);
    if ((OutHashEntry = (this->*GetHashEntryCart)(x, y, theta, cartangle)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntryCart)(x, y, theta, cartangle);
    }

    // need to recompute start heuristics?
    if (EnvNAVXYTHETALATCART.startstateid != OutHashEntry->stateID) {
        bNeedtoRecomputeStartHeuristics = true;
        // because termination condition can be not all states TODO - make it dependent on term. condition
        bNeedtoRecomputeGoalHeuristics = true;
    }

    // set start
    EnvNAVXYTHETALATCART.startstateid = OutHashEntry->stateID;
    EnvNAVXYTHETALATCartCfg.StartX_c = x;
    EnvNAVXYTHETALATCartCfg.StartY_c = y;
    EnvNAVXYTHETALATCartCfg.StartTheta = theta;
    EnvNAVXYTHETALATCartCfg.StartCartAngle = cartangle;

    return EnvNAVXYTHETALATCART.startstateid;
}

void EnvironmentNAVXYTHETALAT::PrintState(
    int stateID,
    bool bVerbose,
    FILE* fOut)
{
    SBPL_INFO("StateID2CoordTableCart Size: %d", StateID2CoordTableCart.size());

#if DEBUG
    if (stateID >= (int)StateID2CoordTableCart.size()) {
        SBPL_ERROR("ERROR in EnvNAVXYTHETALAT... function: stateID illegal (2)\n");
        throw SBPL_Exception();
    }
#endif

    if (fOut == NULL) {
        fOut = stdout;
    }

    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[stateID];

    if (stateID == EnvNAVXYTHETALATCART.goalstateid && bVerbose) {
        SBPL_FPRINTF(fOut, "the state is a goal state\n");
    }

    if (bVerbose) {
        SBPL_FPRINTF(fOut, "X=%d Y=%d Theta=%d CartAngle=%d\n", HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle);
    }
    else
    {
        SBPL_FPRINTF(fOut, "%.3f %.3f %.3f %.3f\n", DISCXY2CONT(HashEntry->X, EnvNAVXYTHETALATCartCfg.cellsize_m), DISCXY2CONT(HashEntry->Y, EnvNAVXYTHETALATCartCfg.cellsize_m), DiscTheta2ContNew(HashEntry->Theta), CartDiscTheta2Cont(HashEntry->CartAngle, CART_THETADIRS));
    }
}

void EnvironmentNAVXYTHETALAT::PrintHashTableHist(FILE* fOut)
{
    int s0 = 0, s1 = 0, s50 = 0, s100 = 0, s200 = 0, s300 = 0, slarge = 0;

    for (int j = 0; j < HashTableSizeCart; j++) {
        if ((int)Coord2StateIDHashTableCart[j].size() == 0)
            s0++;
        else if ((int)Coord2StateIDHashTableCart[j].size() < 5)
            s1++;
        else if ((int)Coord2StateIDHashTableCart[j].size() < 25)
            s50++;
        else if ((int)Coord2StateIDHashTableCart[j].size() < 50)
            s100++;
        else if ((int)Coord2StateIDHashTableCart[j].size() < 100)
            s200++;
        else if ((int)Coord2StateIDHashTableCart[j].size() < 400)
            s300++;
        else
            slarge++;
    }
    SBPL_FPRINTF(fOut, "hash table histogram: 0:%d, <5:%d, <25:%d, <50:%d, <100:%d, <400:%d, >400:%d\n", s0, s1, s50,
                 s100, s200, s300, slarge);
}

int EnvironmentNAVXYTHETALAT::GetFromToHeuristic(int FromStateID, int ToStateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (FromStateID >= (int)StateID2CoordTableCart.size() ||
        ToStateID >= (int)StateID2CoordTableCart.size())
    {
        SBPL_ERROR("ERROR in EnvNAVXYTHETALAT... function: stateID illegal\n");
        throw SBPL_Exception();
    }
#endif

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* FromHashEntry = StateID2CoordTableCart[FromStateID];
    EnvNAVXYTHETALATCARTHashEntry_t* ToHashEntry = StateID2CoordTableCart[ToStateID];

    // TODO - check if one of the gridsearches already computed and then use it.

    return (int)(NAVXYTHETALAT_COSTMULT_MTOMM *
            EuclideanDistance_m(FromHashEntry->X, FromHashEntry->Y, ToHashEntry->X, ToHashEntry->Y) /
            EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs);

}

int EnvironmentNAVXYTHETALAT::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTableCart.size()) {
        throw SBPL_Exception("ERROR in EnvNAVXYTHETALAT... function: stateID illegal");
    }
#endif

    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[stateID];
    // computes distances from start state that is grid2D, so it is EndX_c EndY_c
    int h2D = grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(NAVXYTHETALAT_COSTMULT_MTOMM *
            EuclideanDistance_m(HashEntry->X, HashEntry->Y, EnvNAVXYTHETALATCartCfg.EndX_c, EnvNAVXYTHETALATCartCfg.EndY_c));

    // define this function if it is used in the planner (heuristic backward search would use it)
    return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs);
}

int EnvironmentNAVXYTHETALAT::GetStartHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTableCart.size()) {
        throw SBPL_Exception("ERROR in EnvNAVXYTHETALAT... function: stateID illegal");
    }
#endif

    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[stateID];
    int h2D = grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(NAVXYTHETALAT_COSTMULT_MTOMM *
            EuclideanDistance_m(EnvNAVXYTHETALATCartCfg.StartX_c, EnvNAVXYTHETALATCartCfg.StartY_c, HashEntry->X, HashEntry->Y));

    // define this function if it is used in the planner (heuristic backward search would use it)
    return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETALATCartCfg.nominalvel_mpersecs);
}

int EnvironmentNAVXYTHETALAT::SizeofCreatedEnv()
{
    return (int)StateID2CoordTableCart.size();
}

const EnvNAVXYTHETALATCARTHashEntry_t*
EnvironmentNAVXYTHETALAT::GetStateEntry(int state_id) const
{
    if (state_id >= 0 && state_id < (int)StateID2CoordTableCart.size()) {
        return StateID2CoordTableCart[state_id];
    }
    else {
        return NULL;
    }
}
// end of comments

// for trailer
EnvNAVXYTHETALATCARTHashEntry_t* EnvironmentNAVXYTHETALAT::GetHashEntryCart_lookup(
    int X, int Y, int Theta, int CartAngle)
{
    if (X < 0 || X >= EnvNAVXYTHETALATCartCfg.EnvWidth_c ||
        Y < 0 || Y >= EnvNAVXYTHETALATCartCfg.EnvHeight_c ||
        Theta < 0 || Theta >= EnvNAVXYTHETALATCartCfg.NumThetaDirs)
    {
        return NULL;
    }
    int index = XYTHETACART2INDEX(X,Y,Theta, CartAngle);
    return Coord2StateIDHashTableCart_lookup[index];
}

EnvNAVXYTHETALATCARTHashEntry_t*
EnvironmentNAVXYTHETALAT::GetHashEntryCart_hash(int X, int Y, int Theta, int CartAngle)
{
#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    int binid = GETHASHBIN(X, Y, Theta, CartAngle);
    // SBPL_PRINTF("The size of the binV is %d", (int)Coord2StateIDHashTableCart[binid].size());

// #if DEBUG
//     if ((int)Coord2StateIDHashTableCart[binid].size() > 5) {
//         SBPL_FPRINTF(fDeb, "WARNING: Hash table has a bin %d (X=%d Y=%d) of size %d\n", 
//             binid, X, Y, (int)Coord2StateIDHashTableCart[binid].size());

//         PrintHashTableHist(fDeb);
//     }
// #endif

    // iterate over the states in the bin and select the perfect match
    std::vector<EnvNAVXYTHETALATCARTHashEntry_t*>* binV = &Coord2StateIDHashTableCart[binid];
    for (int ind = 0; ind < (int)binV->size(); ind++) {
        EnvNAVXYTHETALATCARTHashEntry_t* hashentry = binV->at(ind);

        if (hashentry->X == X && hashentry->Y == Y && hashentry->Theta == Theta && hashentry->CartAngle == CartAngle) {
#if TIME_DEBUG
            time_gethash += clock()-currenttime;
#endif
            return hashentry;
        }
    }

#if TIME_DEBUG
    time_gethash += clock()-currenttime;
#endif

    return NULL;
}

EnvNAVXYTHETALATCARTHashEntry_t*
EnvironmentNAVXYTHETALAT::CreateNewHashEntryCart_lookup(int X, int Y, int Theta, int CartAngle)
{
    int i;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = new EnvNAVXYTHETALATCARTHashEntry_t;

    HashEntry->X = X;
    HashEntry->Y = Y;
    HashEntry->Theta = Theta;
    HashEntry->CartAngle = CartAngle;
    HashEntry->iteration = 0;

    HashEntry->stateID = StateID2CoordTableCart.size();

    // insert into the tables
    StateID2CoordTableCart.push_back(HashEntry);

    int index = XYTHETACART2INDEX(X,Y,Theta,CartAngle);

#if DEBUG
    if (Coord2StateIDHashTableCart_lookup[index] != NULL) {
        throw SBPL_Exception("ERROR: creating hash entry for non-NULL hashentry");
    }
#endif

    Coord2StateIDHashTableCart_lookup[index] = HashEntry;

    // insert into and initialize the mappings
    int* entry = new int[NUMOFINDICES_STATEID2IND];
    StateID2IndexMapping.push_back(entry);
    for (i = 0; i < NUMOFINDICES_STATEID2IND; i++) {
        StateID2IndexMapping[HashEntry->stateID][i] = -1;
    }

    if (HashEntry->stateID != (int)StateID2IndexMapping.size() - 1) {
        throw SBPL_Exception("ERROR in Env... function: last state has incorrect stateID");
    }

#if TIME_DEBUG
    time_createhash += clock()-currenttime;
#endif

    return HashEntry;
}

EnvNAVXYTHETALATCARTHashEntry_t*
EnvironmentNAVXYTHETALAT::CreateNewHashEntryCart_hash(int X, int Y, int Theta, int CartAngle)
{
    int i;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = new EnvNAVXYTHETALATCARTHashEntry_t;

    HashEntry->X = X;
    HashEntry->Y = Y;
    HashEntry->Theta = Theta;
    HashEntry->CartAngle = CartAngle;
    HashEntry->iteration = 0;

    HashEntry->stateID = StateID2CoordTableCart.size();

    // insert into the tables
    StateID2CoordTableCart.push_back(HashEntry);

    // get the hash table bin
    i = GETHASHBIN(HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle);

    // insert the entry into the bin
    Coord2StateIDHashTableCart[i].push_back(HashEntry);

    // insert into and initialize the mappings
    int* entry = new int[NUMOFINDICES_STATEID2IND];
    StateID2IndexMapping.push_back(entry);
    for (i = 0; i < NUMOFINDICES_STATEID2IND; i++) {
        StateID2IndexMapping[HashEntry->stateID][i] = -1;
    }

    if (HashEntry->stateID != (int)StateID2IndexMapping.size() - 1) {
        throw SBPL_Exception("ERROR in Env... function: last state has incorrect stateID");
    }

#if TIME_DEBUG
    time_createhash += clock() - currenttime;
#endif

    return HashEntry;
}

// changed for trailer configuration
void EnvironmentNAVXYTHETALAT::GetSuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionV)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // clear the successor array
    SuccIDV->clear();
    CostV->clear();
    SuccIDV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    CostV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    if (actionV != NULL) {
        actionV->clear();
        actionV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    }

    // goal state should be absorbing
    if (SourceStateID == EnvNAVXYTHETALATCART.goalstateid) return;

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[SourceStateID];

    // iterate through actions
    for (aind = 0; aind < EnvNAVXYTHETALATCartCfg.actionwidth; aind++) {
        EnvNAVXYTHETALATCARTAction_t* nav3daction = &EnvNAVXYTHETALATCartCfg.ActionsV[(unsigned int)HashEntry->Theta][(unsigned int)HashEntry->CartAngle][aind];
        // SBPL_INFO("actionwidth %d, aind%d, (%d, %d, %d, %d), (%d, %d, %d, %d)",EnvNAVXYTHETALATCartCfg.actionwidth, aind, HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle, nav3daction->dX, nav3daction->dY, nav3daction->endtheta, nav3daction->endcartangle);
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);
        // TODO
        int newCartAngle = normalizeDiscAngle(nav3daction->endcartangle);

        // SBPL_INFO("new Value(%d, %d, %d, %d)",newX, newY, newTheta, newCartAngle);
        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // get cost
        int cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle, nav3daction);
        // SBPL_INFO("Cost(%d)",cost);
        if (cost >= INFINITECOST) {
            continue;
        }

        EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntryCart)(newX, newY, newTheta, newCartAngle)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntryCart)(newX, newY, newTheta, newCartAngle);
        }

        SuccIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
        if (actionV != NULL) {
            actionV->push_back(nav3daction);
        }
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

// Changed for trailer configuration
void EnvironmentNAVXYTHETALAT::GetPreds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV)
{
    //TODO- to support tolerance, need:
    // a) generate preds for goal state based on all possible goal state variable settings,
    // b) change goal check condition in gethashentry c) change
    //    getpredsofchangedcells and getsuccsofchangedcells functions

    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[TargetStateID];

    // clear the successor array
    PredIDV->clear();
    CostV->clear();
    PredIDV->reserve(EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());
    CostV->reserve(EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());

    // iterate through actions
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionsV = &EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta];
    for (aind = 0; aind < (int)EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size(); aind++) {

        EnvNAVXYTHETALATCARTAction_t* nav3daction = actionsV->at(aind);

        int predX = HashEntry->X - nav3daction->dX;
        int predY = HashEntry->Y - nav3daction->dY;
        int predTheta = nav3daction->starttheta;
        int predCartAngle = nav3daction->startcartangle;


        // skip the invalid cells
        if (!IsValidCell(predX, predY)) {
            continue;
        }

        // get cost
        int cost = GetActionCost(predX, predY, predTheta, predCartAngle, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntryCart)(predX, predY, predTheta, predCartAngle)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntryCart)(predX, predY, predTheta, predCartAngle);
        }

        PredIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

// changed for trailer configuration
void EnvironmentNAVXYTHETALAT::SetAllActionsandAllOutcomes(CMDPSTATE* state)
{
    int cost;

#if DEBUG
    if (state->StateID >= (int)StateID2CoordTableCart.size()) {
        throw SBPL_Exception("ERROR in Env... function: stateID illegal");
    }

    if ((int)state->Actions.size() != 0) {
        throw SBPL_Exception("ERROR in Env_setAllActionsandAllOutcomes: actions already exist for the state");
    }
#endif

    // goal state should be absorbing
    if (state->StateID == EnvNAVXYTHETALATCART.goalstateid) {
        return;
    }

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[state->StateID];

    // iterate through actions
    for (int aind = 0; aind < EnvNAVXYTHETALATCartCfg.actionwidth; aind++) {
        EnvNAVXYTHETALATCARTAction_t* nav3daction = &EnvNAVXYTHETALATCartCfg.ActionsV[(unsigned int)HashEntry->Theta][(unsigned int)HashEntry->CartAngle][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);
        int newCartAngle = NORMALIZEDISCTHETA(nav3daction->endcartangle, CART_THETADIRS);	

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // get cost
        cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        // add the action
        CMDPACTION* action = state->AddAction(aind);

#if TIME_DEBUG
        clock_t currenttime = clock();
#endif

        EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntryCart)(newX, newY, newTheta, newCartAngle)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntryCart)(newX, newY, newTheta, newCartAngle);
        }
        action->AddOutcome(OutHashEntry->stateID, cost, 1.0);

#if TIME_DEBUG
        time3_addallout += clock()-currenttime;
#endif
    }
}

void EnvironmentNAVXYTHETALAT::GetPredsofChangedEdges(
    std::vector<nav2dcell_t> const* changedcellsV,
    std::vector<int>* preds_of_changededgesIDV)
{
    nav2dcell_t cell;
	EnvNAVXYTHETALATCART3Dcell_t affectedcell;
	EnvNAVXYTHETALATCARTHashEntry_t* affectedHashEntry;
  
	//increment iteration for processing savings
	iteration++;
  
	for(int i = 0; i < (int)changedcellsV->size(); i++) 
	{
		cell = changedcellsV->at(i);
		
		//now iterate over all states that could potentially be affected
		for(int sind = 0; sind < (int)affectedpredstatesVc.size(); sind++)
		{
			affectedcell = affectedpredstatesVc.at(sind);

			//translate to correct for the offset
			affectedcell.x = affectedcell.x + cell.x;
			affectedcell.y = affectedcell.y + cell.y;

			//insert only if it was actually generated
      affectedHashEntry = (this->*GetHashEntryCart)(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.cartangle);
			if(affectedHashEntry != NULL && affectedHashEntry->iteration < iteration)
			{
				preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}

    // nav2dcell_t cell;
    // sbpl_xy_theta_cell_t affectedcell;
    // EnvNAVXYTHETALATHashEntry_t* affectedHashEntry;

    // // increment iteration for processing savings
    // iteration++;

    // for (int i = 0; i < (int)changedcellsV->size(); i++) {
    //     cell = changedcellsV->at(i);

    //     // now iterate over all states that could potentially be affected
    //     for (int sind = 0; sind < (int)affectedpredstatesV.size(); sind++) {
    //         affectedcell = affectedpredstatesV.at(sind);

    //         // translate to correct for the offset
    //         affectedcell.x = affectedcell.x + cell.x;
    //         affectedcell.y = affectedcell.y + cell.y;

    //         // insert only if it was actually generated
    //         affectedHashEntry = (this->*GetHashEntry)(affectedcell.x, affectedcell.y, affectedcell.theta);
    //         if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
    //             preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
    //             affectedHashEntry->iteration = iteration; // mark as already inserted
    //         }
    //     }
    // }
}

void EnvironmentNAVXYTHETALAT::GetSuccsofChangedEdges(
    std::vector<nav2dcell_t> const* changedcellsV,
    std::vector<int>* succs_of_changededgesIDV)
{
    nav2dcell_t cell;
	EnvNAVXYTHETALATCART3Dcell_t affectedcell;
	EnvNAVXYTHETALATCARTHashEntry_t* affectedHashEntry;

	throw SBPL_Exception("ERROR: getsuccs is not supported currently");

	//increment iteration for processing savings
	iteration++;

	//TODO - check
	for(int i = 0; i < (int)changedcellsV->size(); i++) 
	{
		cell = changedcellsV->at(i);
		
		//now iterate over all states that could potentially be affected
		for(int sind = 0; sind < (int)affectedsuccstatesVc.size(); sind++)
		{
			affectedcell = affectedsuccstatesVc.at(sind);
      
			//translate to correct for the offset
			affectedcell.x = affectedcell.x + cell.x;
			affectedcell.y = affectedcell.y + cell.y;
      
			//insert only if it was actually generated
      affectedHashEntry = (this->*GetHashEntryCart)(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.cartangle);
			if(affectedHashEntry != NULL && affectedHashEntry->iteration < iteration)
			{
				succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}


    // nav2dcell_t cell;
    // sbpl_xy_theta_cell_t affectedcell;
    // EnvNAVXYTHETALATHashEntry_t* affectedHashEntry;

    // throw SBPL_Exception("ERROR: getsuccs is not supported currently");

    // // increment iteration for processing savings
    // iteration++;

    // // TODO - check
    // for (int i = 0; i < (int)changedcellsV->size(); i++) {
    //     cell = changedcellsV->at(i);

    //     // now iterate over all states that could potentially be affected
    //     for (int sind = 0; sind < (int)affectedsuccstatesV.size(); sind++) {
    //         affectedcell = affectedsuccstatesV.at(sind);

    //         // translate to correct for the offset
    //         affectedcell.x = affectedcell.x + cell.x;
    //         affectedcell.y = affectedcell.y + cell.y;

    //         // insert only if it was actually generated
    //         affectedHashEntry = (this->*GetHashEntry)(affectedcell.x, affectedcell.y, affectedcell.theta);
    //         if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
    //             succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
    //             // mark as already inserted
    //             affectedHashEntry->iteration = iteration;
    //         }
    //     }
    // }
}

void EnvironmentNAVXYTHETALAT::InitializeEnvironmentCart()
{
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry;

    int maxsize = EnvNAVXYTHETALATCartCfg.EnvWidth_c * EnvNAVXYTHETALATCartCfg.EnvHeight_c * EnvNAVXYTHETALATCartCfg.NumThetaDirs*CART_THETADIRS;//Needs to be changed
    SBPL_PRINTF("The parameter = %d, %d, %d. maxsize is %d, %d\n",
         EnvNAVXYTHETALATCartCfg.EnvWidth_c, EnvNAVXYTHETALATCartCfg.EnvHeight_c, EnvNAVXYTHETALATCartCfg.NumThetaDirs,
         maxsize, maxsize <= SBPL_XYTHETALAT_MAXSTATESFORLOOKUP);


    if (maxsize <= SBPL_XYTHETALAT_MAXSTATESFORLOOKUP) {
        SBPL_PRINTF("environment stores states in lookup table\n");

        Coord2StateIDHashTableCart_lookup = new EnvNAVXYTHETALATCARTHashEntry_t*[maxsize];
        for (int i = 0; i < maxsize; i++) {
            Coord2StateIDHashTableCart_lookup[i] = NULL;
        }
        GetHashEntryCart = &EnvironmentNAVXYTHETALAT::GetHashEntryCart_lookup;
        CreateNewHashEntryCart = &EnvironmentNAVXYTHETALAT::CreateNewHashEntryCart_lookup;

        // not using hash table
        HashTableSizeCart = 0;
        Coord2StateIDHashTableCart = NULL;
    }
    else {
        SBPL_PRINTF("environment stores states in hashtable\n");

        // initialize the map from Coord to StateID
        HashTableSizeCart = 4 * 1024 * 1024; // should be power of two
        Coord2StateIDHashTableCart = new std::vector<EnvNAVXYTHETALATCARTHashEntry_t*>[HashTableSizeCart];
        GetHashEntryCart = &EnvironmentNAVXYTHETALAT::GetHashEntryCart_hash;
        CreateNewHashEntryCart = &EnvironmentNAVXYTHETALAT::CreateNewHashEntryCart_hash;

        // not using hash
        Coord2StateIDHashTableCart_lookup = NULL;
    }

    // initialize the map from StateID to Coord
    StateID2CoordTableCart.clear();

    // create start state
    if (NULL == (HashEntry = (this->*GetHashEntryCart)(
            EnvNAVXYTHETALATCartCfg.StartX_c,
            EnvNAVXYTHETALATCartCfg.StartY_c,
            EnvNAVXYTHETALATCartCfg.StartTheta,
            EnvNAVXYTHETALATCartCfg.StartCartAngle)))
    {
        // have to create a new entry
        HashEntry = (this->*CreateNewHashEntryCart)(
                EnvNAVXYTHETALATCartCfg.StartX_c,
                EnvNAVXYTHETALATCartCfg.StartY_c,
                EnvNAVXYTHETALATCartCfg.StartTheta,
                EnvNAVXYTHETALATCartCfg.StartCartAngle);
    }
    EnvNAVXYTHETALATCART.startstateid = HashEntry->stateID;

    // create goal state
    if ((HashEntry = (this->*GetHashEntryCart)(
            EnvNAVXYTHETALATCartCfg.EndX_c,
            EnvNAVXYTHETALATCartCfg.EndY_c,
            EnvNAVXYTHETALATCartCfg.EndTheta,
            EnvNAVXYTHETALATCartCfg.EndCartAngle)) == NULL)
    {
        // have to create a new entry
        HashEntry = (this->*CreateNewHashEntryCart)(
                EnvNAVXYTHETALATCartCfg.EndX_c,
                EnvNAVXYTHETALATCartCfg.EndY_c,
                EnvNAVXYTHETALATCartCfg.EndTheta,
                EnvNAVXYTHETALATCartCfg.EndCartAngle);
    }
    EnvNAVXYTHETALATCART.goalstateid = HashEntry->stateID;

    // initialized
    EnvNAVXYTHETALATCART.bInitialized = true;
}

// examples of hash functions: map state coordinates onto a hash value
// #define GETHASHBIN(X, Y) (Y*WIDTH_Y+X)
// here we have state coord: <X1, X2, X3, X4>
// unsigned int EnvironmentNAVXYTHETALAT::GETHASHBIN(
//     unsigned int X1,
//     unsigned int X2,
//     unsigned int Theta)
// {
//     // SBPL_PRINTF("x, y, theta, cartangle: %d, %d, %d. The result is %d.\n", X1, X2, Theta, 
//     //     inthash(inthash(X1)+(inthash(X2)<<1)+(inthash(Theta)<<2)) & (HashTableSize-1));
//     // SBPL_PRINTF("The value of the inthash: %d, %d, %d, %d, %d\n",
//     //      inthash(X1), inthash(X2), (inthash(X2)<<1), inthash(Theta), (inthash(Theta)<<2));
//     return inthash(inthash(X1) + (inthash(X2) << 1) + (inthash(Theta) << 2)) & (HashTableSize - 1);
// }

unsigned int EnvironmentNAVXYTHETALAT::GETHASHBIN(
    unsigned int X1, 
    unsigned int X2, 
    unsigned int Theta, 
    unsigned int CartAngle)
{
    // SBPL_PRINTF("x, y, theta, cartangle: %d, %d, %d, %d. The result is %d.\n", X1, X2, Theta, CartAngle, 
    //     inthash(inthash(X1)+(inthash(X2)<<1)+(inthash(Theta)<<2)+(inthash(CartAngle)<<3)) & (HashTableSizeCart-1));
    // SBPL_PRINTF("The value of the inthash: %d, %d, %d, %d, %d\n",
    //      inthash(X1), inthash(X2), (inthash(X2)<<1), inthash(Theta), (inthash(Theta)<<2));
    return inthash(inthash(X1)+(inthash(X2)<<1)+(inthash(Theta)<<2)+(inthash(CartAngle)<<3)) & (HashTableSizeCart-1);

	// return inthash(inthash(X1)+(inthash(X2)<<1)+(inthash(Theta)<<2)+(inthash(CartAngle)<<3)) & (HashTableSize-1);
}



//------------------------------------------------------------------------------


void EnvironmentNAVXYTHETALAT::GetLazySuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost,
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionV)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // clear the successor array
    SuccIDV->clear();
    CostV->clear();
    SuccIDV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    CostV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    isTrueCost->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    if (actionV != NULL) {
        actionV->clear();
        actionV->reserve(EnvNAVXYTHETALATCartCfg.actionwidth);
    }

    // goal state should be absorbing
    if (SourceStateID == EnvNAVXYTHETALATCART.goalstateid) {
        return;
    }

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[SourceStateID];

    // iterate through actions
    for (aind = 0; aind < EnvNAVXYTHETALATCartCfg.actionwidth; aind++) {
        EnvNAVXYTHETALATCARTAction_t* nav3daction = &EnvNAVXYTHETALATCartCfg.ActionsV[(unsigned int)HashEntry->Theta][(unsigned int)HashEntry->CartAngle][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);
        int newCartAngle = normalizeDiscAngle(nav3daction->endcartangle);

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // if we are supposed to return the action, then don't do lazy
        if (!actionV) {
            EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
            if ((OutHashEntry = (this->*GetHashEntryCart)(newX, newY, newTheta, newCartAngle)) == NULL) {
                OutHashEntry = (this->*CreateNewHashEntryCart)(newX, newY, newTheta, newCartAngle);
            }
            SuccIDV->push_back(OutHashEntry->stateID);
            CostV->push_back(nav3daction->cost);
            isTrueCost->push_back(false);
            continue;
        }

        // get cost
        int cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, HashEntry->CartAngle, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntryCart)(newX, newY, newTheta, newCartAngle)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntryCart)(newX, newY, newTheta, newCartAngle);
        }

        SuccIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
        isTrueCost->push_back(true);
        if (actionV != NULL) {
            actionV->push_back(nav3daction);
        }
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

// int EnvironmentNAVXYTHETALAT::GetTrueCost(int parentID, int childID)
// {
//     EnvNAVXYTHETALATHashEntry_t* fromHash = StateID2CoordTable[parentID];
//     EnvNAVXYTHETALATHashEntry_t* toHash = StateID2CoordTable[childID];

//     for (int i = 0; i < EnvNAVXYTHETALATCfg.actionwidth; i++) {
//         EnvNAVXYTHETALATAction_t* nav3daction = &EnvNAVXYTHETALATCfg.ActionsV[(unsigned int)fromHash->Theta][i];
//         int newX = fromHash->X + nav3daction->dX;
//         int newY = fromHash->Y + nav3daction->dY;
//         int newTheta = normalizeDiscAngle(nav3daction->endtheta);

//         // skip the invalid cells
//         if (!IsValidCell(newX, newY)) {
//             continue;
//         }

//         EnvNAVXYTHETALATHashEntry_t* hash;
//         if ((hash = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
//             continue;
//         }
//         if (hash->stateID != toHash->stateID) {
//             continue;
//         }

//         // get cost
//         int cost = GetActionCost(fromHash->X, fromHash->Y, fromHash->Theta, nav3daction);

//         if (cost >= INFINITECOST) {
//             return -1;
//         }
//         return cost;
//     }
//     throw SBPL_Exception("This should never happen! we didn't find the state we need to get the true cost for!");
//     return -1;
// }

void EnvironmentNAVXYTHETALAT::GetSuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionV)
{
    GetSuccs(SourceStateID, SuccIDV, CostV, actionV);
}

void EnvironmentNAVXYTHETALAT::GetLazySuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost,
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionV)
{
    GetLazySuccs(SourceStateID, SuccIDV, CostV, isTrueCost, actionV);
}

bool EnvironmentNAVXYTHETALAT::isGoal(int id)
{
    return EnvNAVXYTHETALATCART.goalstateid == id;
}

void EnvironmentNAVXYTHETALAT::GetLazyPreds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // get X, Y for the state
    EnvNAVXYTHETALATCARTHashEntry_t* HashEntry = StateID2CoordTableCart[TargetStateID];

    // clear the successor array
    PredIDV->clear();
    CostV->clear();
    PredIDV->reserve(EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());
    CostV->reserve(EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());

    // iterate through actions
    std::vector<EnvNAVXYTHETALATCARTAction_t*>* actionsV = &EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta];
    for (aind = 0; aind < (int)EnvNAVXYTHETALATCartCfg.PredActionsV[(unsigned int)HashEntry->Theta].size(); aind++)
    {
        EnvNAVXYTHETALATCARTAction_t* nav3daction = actionsV->at(aind);

        int predX = HashEntry->X - nav3daction->dX;
        int predY = HashEntry->Y - nav3daction->dY;
        int predTheta = nav3daction->starttheta;
        int predCartAngle = normalizeDiscAngle(nav3daction->endcartangle);

        //skip the invalid cells
        if (!IsValidCell(predX, predY)) {
            continue;
        }

        EnvNAVXYTHETALATCARTHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntryCart)(predX, predY, predTheta, predCartAngle)) == NULL) {
            OutHashEntry = (this->*CreateNewHashEntryCart)(predX, predY, predTheta, predCartAngle);
        }

        PredIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(nav3daction->cost);
        isTrueCost->push_back(false);
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentNAVXYTHETALAT::GetPredsWithUniqueIds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV)
{
    GetPreds(TargetStateID, PredIDV, CostV);
}

void EnvironmentNAVXYTHETALAT::GetLazyPredsWithUniqueIds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazyPreds(TargetStateID, PredIDV, CostV, isTrueCost);
}

