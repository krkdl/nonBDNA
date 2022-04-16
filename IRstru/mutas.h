
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

//#define RANDOM_TRACE

#define FOLDER_HUMAN   "/Users/nonBDNA/indata/"
#define FOLDER_MUT_DATA "/Users/nonBDNA/mut_set/"
#define FOLDER_IN_DATA "/Users/nonBDNA/cancers/"

#if defined ( RUN_WITH_REPTIME_SEGMENT )
#define FOLDER_OUT_DATA "/Users/nonBDNA/IRstudy/outd30/"
#else
#define FOLDER_OUT_DATA "/Users/nonBDNA/IRstudy/outd30/"
#endif

#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24
#define NO_RANDOM_CYCLES 10000
//#define NO_INT_ZONESIZE 6       // распределение размеров зон мутаций

//#define NO_CANCER_ID 6
#define NO_CANCER_ID 47
#define NO_ZONES_ID 7
#define NO_RT_SEGMENTS 7

int LoadDNKset( const char *pFolder, const char *fSetName );
int LoadXroZones( const char *pFolder, const char *fSetID );
int loadMutation( const char *CancerID ); 

void makeTCXbody( );
void makeRTSEGment( int noSeg );
int loadRepTimeSet( int indCancID);
void markAllRTtag( );

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );


#endif
