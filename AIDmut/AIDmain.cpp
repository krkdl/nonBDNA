// dnk.cpp : Defines the entry point for the console application.
//

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


MUT_SIGNA   signa_C = { 4, -2, { "AACC", "AACT", "AGCC", "AGCT", "TACC", "TACT", "TGCC", "TGCT" } };
MUT_SIGNA   signa_G = { 4, -1, { "GGTT", "AGTT", "GGCT", "AGCT", "GGTA", "AGTA", "GGCA", "AGCA" } };
int curSigna;


long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//                      MUT_signatura       //ZoneID
{
    
		char buffr[1024];
        vector <XROSOMA>::iterator curXRO;
        int indCancer, indZone;

	if ( argc < 2 ) {
		printf ("\nUsage: %s  'MUT_sign'\n", argv[0]);
		return -1;
	}
    
    for ( curSigna=0; curSigna<NO_SIGNA_ID; curSigna++)
        if ( strcmp(argv[1], signa_C.sigset[curSigna] )==0 )
            break;
    if ( curSigna >= NO_SIGNA_ID ) {
        printf ("INV.MUT_sign_ID = '%s'\n", argv[1]);
        return -1;
    }

    sprintf(buffr, "%smutrace.txt", FOLDER_OUT_DATA);
    tsld = fopen(buffr, "w");

    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )  //made 'markTCX_TAG()'-->SET_TCx_TAG(Xtag)
		return -1;
    
    for ( indZone=0; indZone<NO_ZONES_ID; indZone++)  {
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[indZone], 1 ) < 0 )      // made 'SET_MUZTWO_TAG'
            return -1;
    }
    
    for ( indZone=0; indZone<NO_ZONES_ID-1; indZone++)  {     //NO_ZONES_ID - 'cpg'
        clock_t start_Z = clock();
        if ( indZone > 0 )
            clearMUZONE_TAG( );
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[indZone] ) < 0 )      // made 'SET_MUZONE_TAG(Xtag)'
            return -1;
    
        calcNucInZon(  );   //markAllRTtag( );
        for ( indCancer=0; indCancer<NO_CANCER_ID; indCancer++) {
            clock_t start_C = clock();
            if ( ! vecSAMPL.empty() )
                vecSAMPL.clear();
            if ( loadMutation( indCancer ) < 0 )
                return -1;
#ifdef RUN_WITH_NSBI_37
            sprintf(buffr, "%sns%s_%s_%s.txt", FOLDER_OUT_DATA,
                            CancerID[indCancer], ZoneID[indZone], signa_C.sigset[curSigna] );
#else
            sprintf(buffr, "%s%s_%s_%s.txt", FOLDER_OUT_DATA,
                    CancerID[indCancer], ZoneID[indZone], signa_C.sigset[curSigna] );
#endif
            fRezult = fopen(buffr, "w");
            
#ifdef RUN_WITH_NSBI_37
            sprintf(buffr, "%sns%s_%s_%s_A.txt", FOLDER_OUT_DATA,
                            CancerID[indCancer], ZoneID[indZone], signa_C.sigset[curSigna] );
#else
            sprintf(buffr, "%s%s_%s_%s.txt", FOLDER_OUT_DATA,
                    CancerID[indCancer], ZoneID[indZone], signa_C.sigset[curSigna] );
#endif
            fRezAPO = fopen(buffr, "w");
            
            printTableCap( fRezult );
            printTableCap( fRezAPO );

            makeTCXbody( );
            switchToXRO( ); // pointXRO( );     // segment = whole Xromo
        
            for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
                vecSAMPL[nS].calcREALstat( );
//              vecSAMPL[nS].countAPOBECmut( );
                vecSAMPL[nS].correctREAL_Xmut(nS);
                vecSAMPL[nS].makeRndCycle(nS, CancerID[indCancer]);
                vecSAMPL[nS].printSamplLine( fRezult, 0);
                vecSAMPL[nS].printSamplLine( fRezAPO, 1);
            }
            clock_t finish_C = clock();
            double duration = (double)(finish_C - start_C) / CLOCKS_PER_SEC;
            printf("------Cancer='%s' : dT=%5.2f sec\n", CancerID[indCancer], duration);
            fclose (fRezult);
            fclose (fRezAPO);
        }
        clock_t finish_Z = clock();
        double duration = (double)(finish_Z - start_Z) / CLOCKS_PER_SEC;
        printf("\n==========Zone='%s' : dT=%5.2f sec=============\n", ZoneID[indZone], duration);
    }
    
	return 0;
}
////////////////////////////////////////////////////////

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

