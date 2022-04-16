
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

#define DEBUG_TRACE_TSLD

//#define RUN_WITH_REPTIME_SEGMENT
//#define DEBUG_TRACE_RTSEG

#if defined ( DEBUG_TRACE_RTSEG )
#define DEBUG_TRACE_TSLD
#endif

#define RTS_MODE
#define CALC_OUT_ZONES

#define FOLDER_HUMAN   "/Users/nonBDNA/indata/"
#define FOLDER_IN_DATA "/Users/nonBDNA/cancers/"

#define FOLDER_OUT_DATA "/Users/nonBDNA/GQstudy/RTS/OutZons/"
#define RTS_DATA "/Users/nonBDNA/GQstudy/ReplicationTimingStrand.txt"
#define CPG_ZoneName  "cpgIslandExt.hg19.txt"

#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24

#define NO_CANCER_ID 6
#define NO_ZONES_ID 8
#define NO_RT_SEGMENTS 7

int LoadDNKset( const char *pFolder, const char *fSetName );
int LoadXroZones( const char *pFolder, const char *fSetID, int TagOnly=0 );

int LoadCPGzone( const char *pFolder, const char *cpgZoneName, int TagOnly=0 );
int loadMutation( int indCancer );
int readRTSzones( const char *fRTSid );
void makeRTSEGment( int noSeg );

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );


#endif
