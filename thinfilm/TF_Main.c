/***************************************************************************
 *
 *  Function    : TF Main
 *  Description : main routine for ParaDiS/ thinfilm simulation
 *
 **************************************************************************/
#include <stdio.h>
#include <time.h>
#include "Home.h"

#ifdef _THINFILM
#include "TF.h"
#endif

#include "Init.h"


#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef FPES_ON
#include <fpcontrol.h>
#endif

void ParadisInit(int argc, char *argv[],Home_t **homeptr);

#ifdef _THINFILM
#include "TF.h"
void ParadisStep(Home_t *home,ThinFilm_t *thinfilm);
#else
void ParadisStep(Home_t *home);
#endif



void ParadisFinish(Home_t *home);

main (int argc, char *argv[])
{
        int     cycleEnd, memSize, initialDLBCycles;
        time_t  tp;
        Home_t  *home;
        Param_t *param;


#ifdef _THINFILM
	ThinFilm_t *thinfilm;
#endif
/*
 *      On some systems, the getrusage() call made by Meminfo() to get
 *      the memory resident set size does not work properly.  In those
 *      cases, the function will try to return the current heap size 
 *      instead.  This initial call allows meminfo() to get a copy of
 *      the original heap pointer so subsequent calls can calculate the
 *      heap size by taking the diference of the original and current
 *      heap pointers.
 */
        Meminfo(&memSize);

/*
 *      on linux systems (e.g. MCR) if built to have floating point exceptions
 *      turned on, invoke macro to do so
 */
   
#ifdef FPES_ON
        unmask_std_fpes();
#endif

        ParadisInit(argc, argv, &home);
   
#ifdef _THINFILM    
        TF_Init(home,&thinfilm);
#endif

        home->cycle      = home->param->cycleStart;

        param            = home->param;
        cycleEnd         = param->cycleStart + param->maxstep;
        initialDLBCycles = param->numDLBCycles;

/*
 *      Perform the needed number (if any) of load-balance-only
 *      steps before doing the main processing loop.  These steps
 *      perform only the minimal amount of stuff needed to
 *      estimate per-process load, move boundaries and migrate
 *      nodes among processsors to get a good initial balance.
 */
        TimerStart(home, INITIALIZE);

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Beginning %d load-balancing steps at %s",
                   initialDLBCycles, asctime(localtime(&tp)));
        }

        while (param->numDLBCycles > 0) {
#ifdef _THINFILM
	  ParadisStep(home,thinfilm);
#else
	  ParadisStep(home);
#endif
            home->cycle++;
            param->numDLBCycles--;
        }

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Completed load-balancing steps at %s",
                   asctime(localtime(&tp)));
        }

        TimerStop(home, INITIALIZE);

/*
 *      Any time spent doing the initial DLB-only steps should
 *      just be attributed to initialization time, so be sure to
 *      reset the other timers before going into the main
 *      computational loop
 */
        TimerInitDLBReset(home);

/*
 *      The cycle number may have been incremented during the initial
 *      load-balance steps, so reset it to the proper starting
 *      value before entering the main processing loop.
 */
        home->cycle = home->param->cycleStart;

/*      For checking only */
#ifdef _THINFILM
#ifdef _PRINTSTRESS
        CellCharge(home);
	SortNativeNodes(home);

#ifdef _TFIMGSTRESS
	TF_Step(home,thinfilm);
#endif
	PrintStress(home,thinfilm);
#endif
#endif



        while (home->cycle < cycleEnd) {

#ifdef _THINFILM
#ifdef _TFIMGSTRESS
	  SortNativeNodes(home);
	  CellCharge(home); 
	  TF_Step(home,thinfilm);
	  CommSendGhosts(home);
#endif
#endif


#ifdef _THINFILM
#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
	  BuildSurfaceSegList(home, thinfilm);
#endif
	  ParadisStep(home,thinfilm);
#else
	  ParadisStep(home);
#endif
	  
          TimerClearAll(home);
        }

        ParadisFinish(home);


#ifdef _THINFILM
#ifdef  _TFIMGSTRESS
        TF_Finish(thinfilm);
	FreeCellCters();
#endif
#endif

        exit(0);
}
