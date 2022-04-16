// dnk.cpp : Defines the entry point for the console application.
//
//#include <iostream>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	FILE *tsld=NULL;
#endif

	FILE *fRezult=NULL;
    FILE *fRezAPO=NULL;

const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z", "cpg"};

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

//void testBVin ();

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//                   "mutation_File_Path" 
{
    
		char buffr[1024];
        vector <XROSOMA>::iterator curXRO;
        int indCancer, indZone;

	if ( argc < 3 ) {
		printf ("\nUsage: %s  'Cancer_ID' 'ZoneID'\n", argv[0]);
		return -1;
	}
    
    for ( indCancer=0; indCancer<NO_CANCER_ID; indCancer++)
        if ( strcmp(argv[1], CancerID[indCancer])==0 )
            break;
    if ( indCancer >=NO_CANCER_ID ) {
        printf ("INV.Cancer_ID = '%s'\n", argv[1]);
        return -1;
    }
    for ( indZone=0; indZone<NO_ZONES_ID; indZone++)
        if ( strcmp(argv[2], ZoneID[indZone])==0 )
            break;
    if ( indZone >=NO_ZONES_ID ) {
        printf ("INV.Zone_ID = '%s'\n", argv[2]);
        return -1;
    }
    
#ifdef DEBUG_TRACE_TSLD
//    sprintf(buffr, "%strace_%s_%s.txt", FOLDER_OUT_DATA,argv[1],argv[2]);
    sprintf(buffr, "%smutrace.txt", FOLDER_OUT_DATA);
    tsld = fopen(buffr, "w");
#endif
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )  //made 'markTCX_TAG()'-->SET_TCx_TAG(Xtag)
		return -1;
    
    if ( loadRepTimeSet(indCancer) < 0 )
        return -1;
    
    for ( int Z=0; Z<NO_ZONES_ID; Z++)  {
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[Z], 1 ) < 0 )      // made 'SET_MUZTWO_TAG'
            return -1;
    }
    if ( LoadXroZones( FOLDER_IN_DATA, argv[2] ) < 0 )      // made 'SET_MUZONE_TAG(Xtag)'
        return -1;
    
	if ( loadMutation( indCancer ) < 0 )
		return -1;

#if defined ( RUN_WITH_REPTIME_SEGMENT )
    clock_t start0 = clock();
    for ( int nSegm=0;  nSegm<NO_RT_SEGMENTS; nSegm++ ) {
        sprintf(buffr, "%s%s_%s_RT%d_t.txt", FOLDER_OUT_DATA, argv[1], argv[2], nSegm);
            fRezult = fopen(buffr, "w");
        sprintf(buffr, "%s%s_%s_RT%d_apo_t.txt", FOLDER_OUT_DATA, argv[1], argv[2], nSegm);
            fRezAPO = fopen(buffr, "w");
        printTableCap( fRezult );
        printTableCap( fRezAPO );
        makeRTSEGment(nSegm);
        makeTCXbody( );
        switchToRT( );   // pointRT( );  //
        
        printf("\nRANDOMcycles (%ld samples) RTseg=%d .....\n", vecSAMPL.size(), nSegm);
        clock_t start1 = clock();
        for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
            vecSAMPL[nS].calcREALstat( );
//            vecSAMPL[nS].countAPOBECmut( );
            vecSAMPL[nS].correctREAL_Xmut(nS);
            vecSAMPL[nS].makeRndCycle(nS, argv[1]);
            vecSAMPL[nS].printSamplLine( fRezult, 0);
            vecSAMPL[nS].printSamplLine( fRezAPO, 1);
        }
        clock_t finish1 = clock();
        double duration = (double)(finish1 - start1) / CLOCKS_PER_SEC;
        printf("====end of SampCycles : dT=%5.2f sec\n", duration);
        fclose (fRezult);
        fclose (fRezAPO);
    }
    clock_t finish0 = clock();
    double duration0 = (double)(finish0 - start0) / CLOCKS_PER_SEC;
    printf("end of SegCycles : dT=%5.2f sec\n", duration0);
#else
    sprintf(buffr, "%s%s_%s.txt", FOLDER_OUT_DATA, argv[1], argv[2]);
    fRezult = fopen(buffr, "w");
    sprintf(buffr, "%s%s_%s_apo.txt", FOLDER_OUT_DATA, argv[1], argv[2] );
    fRezAPO = fopen(buffr, "w");
    printTableCap( fRezult );
    printTableCap( fRezAPO );
    
    markAllRTtag( );
    makeTCXbody( );
    switchToXRO( ); // pointXRO( );     // segment = whole Xromo
    printf("\nRANDOMcycles (%ld samples) .....\n", vecSAMPL.size() );
    clock_t start1 = clock();
    for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
        vecSAMPL[nS].calcREALstat( );
        vecSAMPL[nS].correctREAL_Xmut(nS);
        vecSAMPL[nS].makeRndCycle(nS, argv[1]);
        vecSAMPL[nS].printSamplLine( fRezult, 0);
        vecSAMPL[nS].printSamplLine( fRezAPO, 1);
    }
    clock_t finish1 = clock();
    double duration = (double)(finish1 - start1) / CLOCKS_PER_SEC;
    printf("end of Cycles : dT=%5.2f sec\n", duration);
    fclose (fRezult);
    fclose (fRezAPO);
#endif
    
	return 0;
}
//////////////////////////////////////////////////////////////////////
//====================================================================

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in )
{
    char buff[4096];
    int lng;
    
    if ( fgets(shortRec, sizeRec, f_in)==NULL )
        return 0;
    lng = (int)strlen(shortRec);
    if ( *(shortRec+lng-1) != '\n' )
        while ( 1 ) {
            fgets(buff, sizeof(buff), f_in);
            if ( buff[strlen(buff)-1] == '\n' )
                break;
        }
    
    return lng;
}
////////////////////////////////////////////////////////

