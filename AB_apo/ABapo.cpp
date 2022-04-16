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
const char ZoneID[NO_ZONES_ID][4] = { "apr", "dr", "gq", "ir", "mr", "str", "z"};

#ifdef  ALL_MUTA
char curCancID[8];
#endif

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
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )   // 'markTCX_TAG()
		return -1;

#ifdef  ALL_MUTA
    sprintf(buffr, "%sGQ_AllMu.txt", FOLDER_OUT_DATA);
#else
    sprintf(buffr, "%sGQ_%s.txt", FOLDER_OUT_DATA, argv[1]);
#endif
    
    fRezult = fopen(buffr, "w");
    printTableCap(fRezult );
    
#ifdef  ALL_MUTA
    for ( indCancer=0; indCancer<NO_CANCER_ID; indCancer++)   {
        strcpy(curCancID, CancerID[indCancer]);
        if ( ! vecSAMPL.empty() )
            vecSAMPL.clear();
        if ( loadMutation( indCancer ) < 0 )
            return -1;
        printf("calcAB_APO() Cancer=%s ...\n", CancerID[indCancer] );
        for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
//            printf("calcAB_APO() sampl=%d...\n", nS);
            vecSAMPL[nS].calc_AB_APO(  );
            vecSAMPL[nS].printSampLine( fRezult );
        }
    }
#else
    if ( loadMutation( indCancer ) < 0 )
		return -1;
    printf("calcAB_APO() Cancer=%s ...\n", CancerID[indCancer] );
    for ( int nS=0; nS<vecSAMPL.size(); nS++ ) {
        vecSAMPL[nS].calc_AB_APO(  );
        vecSAMPL[nS].printSampLine( fRezult );
    }
#endif
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

//////////////////////////////////////////////////////////////////////////////////

void printTableCap( FILE *fRez )
{
#ifdef  ALL_MUTA
    fprintf(fRez, "mutID\t#s\tsample\tTrgA\tTrgB\tApoA\tApoB");
#else
    fprintf(fRez, "#s\tsample\tTrgA\tTrgB\tApoA\tApoB");
#endif

    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\tX.%d_TrgA\tX.%d_TrgB\tX.%d_ApoA\tX.%d_ApoB",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");

    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSampLine( FILE *fRez)
{
    
    int sumTrg_A=0, sumTrg_B=0;
    int sumApo_A=0, sumApo_B=0;
    

    for ( int nX=0; nX<vecDNK.size(); nX++ ) {
        sumTrg_A += vecDNK[nX].AcntTrg;
        sumTrg_B += vecDNK[nX].BcntTrg;
        sumApo_A += AcntSamp[nX];
        sumApo_B += BcntSamp[nX];
    }

#ifdef  ALL_MUTA
    fprintf(fRez, "%s\t%d\t%s\t%d\t%d\t%d\t%d",
            curCancID, (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTrg_A, sumTrg_B, sumApo_A, sumApo_B );
#else
    fprintf(fRez, "%d\t%s\t%d\t%d\t%d\t%d",
            (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTrg_A, sumTrg_B, sumApo_A, sumApo_B );
#endif
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%d\t%d\t%d\t%d",
                vecDNK[nX].AcntTrg, vecDNK[nX].BcntTrg,
                AcntSamp[nX], BcntSamp[nX]);
    }
    fprintf(fRez, "\n");

    return;
}
////////////////////////////////////////////////////////
