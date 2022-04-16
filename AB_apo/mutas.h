
#ifndef _MUTAS_H__
#define _MUTAS_H__

#include "vector"
using namespace std;

#define ALL_MUTA

#define DEBUG_TRACE_TSLD


#define FOLDER_HUMAN "/Users/nonBDNA/indata/"
#define FOLDER_IN_DATA "/Users/nonBDNA/cancers/"

#define FOLDER_OUT_DATA "/Users/nonBDNA/AB_apo/outd/"


#define XRO_SET_FILENAME "xromo_set.txt"

#define NO_HUMAN_XRO 24

#define NO_CANCER_ID 6
#define NO_ZONES_ID 7

int LoadDNKset( const char *pFolder, const char *fSetName );
int loadMutation( int indCancer );


int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );


#endif
