#include <sys/types.h>
#include <sys/stat.h>
#include "Home.h"
#include "Util.h"
/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVelocity
 *      Description: Write the velocity of all active nodes to a file.
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *
 *-------------------------------------------------------------------------*/
void WriteVelocity(Home_t *home, char *baseFileName, int ioGroup,
                   int firstInGroup, int writePrologue, int writeEpilogue)
{
        int      i, newNodeKeyPtr,nodecount,nodeid;
        real8    vx, vy, vz,vaver,burgMag;
        char     fileName[256];
        Node_t   *node;
        Param_t  *param = home->param;
        FILE     *fp;
        struct stat statbuf;
		char           workingDir[512];
		(void *)getcwd(workingDir, sizeof(workingDir) - 1);
/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s/%s", workingDir,DIR_VELOCITY,
                     baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s/%s.%d", workingDir,DIR_VELOCITY,
                     baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (firstInGroup) {
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteVelocity: Open error %d on %s\n", errno, fileName);
            }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteVelocity: Open error %d on %s\n", errno, fileName);
            }
        }

/*
 *      Write velocity data for all nodes in this domain
 */
        burgMag = param->burgMag;
        newNodeKeyPtr = home->newNodeKeyPtr;
 		nodecount=0;
 		vaver=0.0;
 		nodeid=-1;
        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            vx =  node->vX * burgMag;
            vy =  node->vY * burgMag;
            vz =  node->vZ * burgMag;
            
			nodecount=nodecount+1;
			vaver=vaver+sqrt(vx*vx+vy*vy+vz*vz);
			
			
			nodeid=node->myTag.domainID;
            fprintf(fp, "%e %e %e %d  (%d,%d, %d)\n", vx, vy, vz, node->sgnv,
                    node->myTag.domainID, node->myTag.index, node->constraint);
        }

        fclose(fp);
         
         //snprintf(fileName, sizeof(fileName), "%s/nopeudet",
                          //DIR_PROPERTIES);
                 //fp = fopen(fileName, "a");
                 //fprintf(fp,"%e %e %d \n", param->timeNow, vaver/nodecount,nodeid);
                 //fclose(fp);
         
         
         
 //printf("TASSA! %e ",vaver/nodecount);
        return;
}
