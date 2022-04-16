//
//  main.cpp
//  irstudy
//
//  Created by Ирина Пономарева on 16/04/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

//#include <iostream>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"
//#include "gar_fu.h"

#ifdef DEBUG_TRACE_TSLD
	FILE *tsld=NULL;
#endif

	FILE *fRezult=NULL;
    FILE *fRezAPO=NULL;

const char CancerID[NO_CANCER_ID][8] = { "BLCA-US", "BOCA-UK", "BRCA-EU", "BRCA-UK", "BRCA-US",
    "BTCA-SG", "CESC-US", "CLLE-ES", "CMDI-UK", "COAD-US",
    "DLBC-US", "EOPC-DE", "ESAD-UK", "GACA-CN", "GBM-US",
    "HNSC-US", "KICH-US", "KIRC-US", "KIRP-US", "LAML-KR",
    "LGG-US",  "LICA-FR", "LIHC-US", "LINC-JP", "LIRI-JP",
    "LUAD-US", "LUSC-US", "MALY-DE", "MELA-AU", "ORCA-IN",
    "OV-AU",   "OV-US",   "PACA-AU", "PACA-CA", "PAEN-AU",
    "PAEN-IT", "PBCA-DE", "PRAD-CA", "PRAD-UK", "PRAD-US",
    "READ-US", "RECA-EU", "SARC-US", "SKCM-US", "STAD-US",
    "THCA-US", "UCEC-US"
};
//const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z"};

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

void printIRCap( FILE *fLoopRez, FILE *fRez );
void printIRsampl( FILE *fLoopRez, FILE *fRez, int nSamp );

void commonTestCnt( );

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
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )  //made 'markTCX_TAG()'-->SET_TCx_TAG(Xtag)
		return -1;

    if ( LoadXroZones( FOLDER_IN_DATA, "ir" ) < 0 )      // made 'SET_MUZONE_TAG' & 'SET_CORN_TAG'
        return -1;
    
	if ( loadMutation( CancerID[indCancer] ) < 0 )    //loadMutation( indCancer )
		return -1;
    
    sprintf(buffr, "%sirstem_%s.txt", FOLDER_OUT_DATA, argv[1]);
    fRezult = fopen(buffr, "w");
    
    sprintf(buffr, "%sirloop_%s.txt", FOLDER_OUT_DATA, argv[1]);
    fRezAPO = fopen(buffr, "w");
    printIRCap( fRezAPO, fRezult );
    
    clock_t start = clock();
    clock_t finish = clock();
    double duration;
    
    for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {    //vecSAMPL.size();
        start = clock();
        printf("clc_IR samp=%d... ", nS);
        
        if ( nS > 0 )
            vecSAMPL[nS-1].markMutations(0);
        vecSAMPL[nS].markMutations(1);
        
        for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
            vecDNK[nX].calcIRstat ( vecSAMPL[nS].vMutNuc[nX] );
        }
        
        printIRsampl( fRezAPO, fRezult, nS );
        
        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        printf("  dT_clc=%5.2fs\n", duration);
    }
    fclose (fRezAPO);
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

