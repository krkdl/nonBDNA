//
//  main.cpp
//  Zstudy
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

//const char CancerID[NO_CANCER_ID][8] = { "BLCA", "BRCA", "CESC", "HNSC", "LUAD", "LUSC"};
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
const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z", "cpg"};

long XROSOMA::maxXsize = 0;
extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//                   "cancerID"
{
    
		char buffr[1024];
//        vector <XROSOMA>::iterator curXRO;
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
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )
		return -1;
    
    for ( int Z=0; Z<NO_ZONES_ID; Z++)  {
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[Z], 1 ) < 0 )      // made 'SET_MUZTWO_TAG'
            return -1;
    }
    if ( LoadXroZones( FOLDER_IN_DATA, "z" ) < 0 )      // made 'SET_MUZONE_TAG(Xtag)'
        return -1;
    
    printf("calcZ_Targets( )  ...\n");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ )
        vecDNK[nX].calcZ_Targets( );
    
	if ( loadMutation( CancerID[indCancer] ) < 0 )
		return -1;
    
    sprintf(buffr, "%sZ_%s.txt", FOLDER_OUT_DATA, argv[1]);
    fRezult = fopen(buffr, "w");
    
    printTableCap(fRezult );
    
    for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
        printf("calcZ_Mutation sampl=%d...\n", nS);
        vecSAMPL[nS].calcZ_Mutations( ); 
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
//////////////////////////////////////////////////////////////////////////////////

void printTableCap( FILE *fRez )
{
    
    fprintf(fRez, "#s\tsample\ttTCin\ttTCout\ttCTin\ttCTout\ttCCin\ttCCout");
    fprintf(fRez,           "\tmTCin\tmTCout\tmCTin\tmCTout\tmCcin\tmCcout\tmcCin\tmcCout");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\tX.%d_tTCin\tX.%d_tTCout\tX.%d_tCTin\tX.%d_tCTout\tX.%d_tCCin\tX.%d_tCCout",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
        fprintf(fRez, "\tX.%d_mTCin\tX.%d_mTCout\tX.%d_mCTin\tX.%d_mCTout\tX.%d_mCcin\tX.%d_mCcout\tX.%d_mcCin\tX.%d_mcCout",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSampLine( FILE *fRez)
{
    CNTs sumTrgIN, sumTrgOUT;
    CNTs sumMutIN, sumMutOUT;
    vector <XROSOMA>::iterator curXr;
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        curXr = vecDNK.begin() + nX;
        sumTrgIN.CT += curXr->cTrgIN.CT;
        sumTrgIN.TC += curXr->cTrgIN.TC;
        sumTrgIN.Cc += curXr->cTrgIN.Cc;
        sumTrgIN.cC += curXr->cTrgIN.cC;
            sumTrgOUT.CT += curXr->cTrgOUT.CT;
            sumTrgOUT.TC += curXr->cTrgOUT.TC;
            sumTrgOUT.Cc += curXr->cTrgOUT.Cc;
            sumTrgOUT.cC += curXr->cTrgOUT.cC;
        sumMutIN.CT += cMutIN[nX].CT;
        sumMutIN.TC += cMutIN[nX].TC;
        sumMutIN.Cc += cMutIN[nX].Cc;
        sumMutIN.cC += cMutIN[nX].cC;
            sumMutOUT.CT += cMutOUT[nX].CT;
            sumMutOUT.TC += cMutOUT[nX].TC;
            sumMutOUT.Cc += cMutOUT[nX].Cc;
            sumMutOUT.cC += cMutOUT[nX].cC;
    }
    
    
    fprintf(fRez, "%d\t%s", (int)(this-&vecSAMPL[0]), SampName.c_str() );
//               tTCin tTCout tCTin tCTout tCCin tCCout mTCin mTCout mCTin mCTout mCcin mCcout mcCin mcCout
    fprintf(fRez, "\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
            sumTrgIN.CT,  sumTrgOUT.CT,  sumTrgIN.TC, sumTrgOUT.TC,
            (sumTrgIN.cC+sumTrgIN.Cc), (sumTrgOUT.cC+sumTrgOUT.Cc),
            sumMutIN.CT, sumMutOUT.CT, sumMutIN.TC, sumMutOUT.TC,
            sumMutIN.Cc, sumMutOUT.Cc, sumMutIN.cC, sumMutOUT.cC );
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                curXr->cTrgIN.CT, curXr->cTrgOUT.CT, curXr->cTrgIN.TC, curXr->cTrgOUT.TC,
                (curXr->cTrgIN.Cc + curXr->cTrgIN.cC), (curXr->cTrgOUT.Cc + curXr->cTrgOUT.cC),
                cMutIN[nX].CT, cMutOUT[nX].CT, cMutIN[nX].TC, cMutOUT[nX].TC,
                cMutIN[nX].Cc, cMutOUT[nX].Cc, cMutIN[nX].cC, cMutOUT[nX].cC );
    }
    fprintf(fRez, "\n");
    
    return;
}
/////////////////////////////////////////////////////////////////////////////////


