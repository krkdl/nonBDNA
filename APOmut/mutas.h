
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

#define RUN_WITH_RANDOM_CYCLES

#define DEBUG_TRACE_TSLD

#define FOLDER_HUMAN "/Users/nonBDNA/indata/"
#define FOLDER_IN_DATA "/Users/nonBDNA/cancers/"
#define FOLDER_MUT_DATA "/Users/nonBDNA/mut_set/"
#define CPG_ZoneName  "cpgIslandExt.hg19.txt"

#if defined ( RUN_WITH_REPTIME_SEGMENT )
#define FOLDER_OUT_DATA "/Users/nonBDNA/APOmut/outRTS/"
#else
#define FOLDER_OUT_DATA "/Users/nonBDNA/APOmut/outd/"
#endif

#if defined ( RUN_WITH_RANDOM_CYCLES )
#define NO_RANDOM_CYCLES 10000
#else
#define NO_RANDOM_CYCLES 0
#endif

#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24

#define NO_CANCER_ID 6
#define NO_ZONES_ID 8
#define NO_RT_SEGMENTS 7

int LoadDNKset( const char *pFolder, const char *fSetName );
int LoadXroZones( const char *pFolder, const char *fSetID, int TagOnly=0 ); //  TagOnly=1 :: SET_MUZTWO_TAG
//  not create 'vecMuZon'
int LoadCPGzone( const char *pFolder, const char *cpgZoneName, int TagOnly=0 );
int loadMutation( int indCancer );

void makeTCXbody( );//selectAllTCx ();
void makeRTSEGment( int noSeg );
int loadRepTimeSet( int indCancID);
void markAllRTtag( );

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );


#endif
