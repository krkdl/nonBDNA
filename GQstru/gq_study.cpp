//
//  main.cpp
//  irstudy
//
//  Created by Ирина Пономарева on 16/04/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
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

const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z", "cpg"}; 

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

void  calcGQ_CrossTest();

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//                   "mutation_File_Path" 
{
    
		char buffr[1024];
        vector <XROSOMA>::iterator curXRO;
        int indCancer;

	if ( argc < 2 ) {
		printf ("\nUsage: %s  'Cancer_ID'\n", argv[0]);
		return -1;
	}
    
    for ( indCancer=0; indCancer<NO_CANCER_ID; indCancer++)
        if ( strcmp(argv[1], CancerID[indCancer])==0 )
            break;
    if ( indCancer >=NO_CANCER_ID ) {
        printf ("INV.Cancer_ID = '%s'\n", argv[1]);
        return -1;
    }
    
#ifdef DEBUG_TRACE_TSLD
    sprintf(buffr, "%smutrace.txt", FOLDER_OUT_DATA);
    tsld = fopen(buffr, "w");
#endif
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )      //  made ' 'markTCX_TAG()
		return -1;
    
#ifdef RTS_MODE
    if ( readRTSzones( RTS_DATA ) <= 0 )
        return -1;
#endif
    
    for ( int Z=0; Z<NO_ZONES_ID; Z++)  {   //NO_ZONES_ID
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[Z], 1 ) < 0 )      // made 'SET_MUZTWO_TAG'
            return -1;
    }
    if ( LoadXroZones( FOLDER_IN_DATA, "gq", 0 ) < 0 )      // made 'SET_MUZONE_TAG' & 'markTCX_TAG()
        return -1;

#ifdef CALC_OUT_ZONES
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ )
            vecDNK[nX].calcTCx_OutZ ( );            //   calc Targets out All Zones
#endif
    
	if ( loadMutation( indCancer ) < 0 )
		return -1;
    

#ifdef CALC_OUT_ZONES
    sprintf(buffr, "%sGQ_%sout.txt", FOLDER_OUT_DATA, argv[1]);
#else
    sprintf(buffr, "%sGQ_%sin.txt", FOLDER_OUT_DATA, argv[1]);
#endif
    fRezult = fopen(buffr, "w");
    printTableCap(fRezult );
    
    clock_t start;
    clock_t finish;
    
    for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
        printf("calc_GQ sampl=%d...\n", nS);
        vecSAMPL[nS].calcGQ_APO(  );
        vecSAMPL[nS].printSampLine( fRezult );
    }
    fclose (fRezult);
    
	return 0;
}
//////////////////////////////////////////////////////////////////////
//====================================================================
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

void calcGQ_CrossTest()
{
    vector <XROSOMA>:: iterator itXR;
    int sumPlus=0, sumMinus=0, sumCross=0;
    fprintf(tsld, "Xr\tStack\tStrand\tCross\n");
    for ( itXR=vecDNK.begin(); itXR!=vecDNK.end(); itXR++ )  {
        int nPlus=0, nMinus=0, nCross=0;
        for (char *pT=itXR->Xtag; pT<(itXR->Xtag+itXR->Xsize); pT++ )      {
            if ( GET_MUZONE_TAG(pT) )   {
                nPlus ++;
                if ( GET_MUZMNS_TAG(pT) )
                    nCross ++;
            }
            else
            if ( GET_MUZMNS_TAG(pT) )
                nMinus ++;
        }
        sumPlus += nPlus;
        sumMinus += nMinus;
        sumCross += nCross;
        fprintf(tsld, "%d\t%d\t%d\t%d\n", itXR->chrNum, nPlus, nMinus, nCross );
    }
    fprintf(tsld, "sum\t%d\t%d\t%d\n", sumPlus, sumMinus, sumCross );
    return;
}
////////////////////////////////////////////////////////


