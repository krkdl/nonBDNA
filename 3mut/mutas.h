
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

//#define RUN_WITH_REPTIME_SEGMENT

#define DEBUG_TRACE_RTSEG

#define DEBUG_TRACE_TSLD

//#define RANDOM_TRACE

#define FOLDER_HUMAN   "/Users/nonBDNA/indata/"
#define FOLDER_IN_DATA "/Users/nonBDNA/cancers/"

#if defined ( RUN_WITH_REPTIME_SEGMENT )
#define FOLDER_OUT_DATA "/Users/nonBDNA/3mut/RTS/outd/"
#else
#define FOLDER_OUT_DATA "/Users/nonBDNA/3mut/outd/"
#endif


#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24
#define NO_RANDOM_CYCLES 10000
//#define NO_INT_ZONESIZE 6       // распределение размеров зон мутаций

#define NO_CANCER_ID 6
#define NO_ZONES_ID 7
#define NO_RT_SEGMENTS 7


int LoadDNKset( const char *pFolder, const char *fSetName );
int LoadXroZones( const char *pFolder, const char *fSetID );
int loadMutation( int indCancer );

void makeTCXbody( );
void makeRTSEGment( int noSeg );
int loadRepTimeSet( int indCancID);
void markAllRTtag( );

//void makeRezultName(const char *f_in, int nS, string &sFiName, const char *zID);
int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );


#endif
