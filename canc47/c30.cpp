//
//  main.cpp
//  cancer30
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

void printIRCap( FILE *fLoopRez, FILE *fRez );
void printIRsampl(FILE *fLoopRez, FILE *fRez, int nSamp ); 

/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) 
//                   "mutation_File_Path" 
{
        char buffr[1024];
        vector <XROSOMA>::iterator curXRO;
        int indZ;

	if ( argc < 2 ) {
		printf ("\nUsage: %s  'Struct_ID'\n", argv[0]);
		return -1;
	}
    
    for ( indZ=0; indZ<NO_ZONES_ID; indZ++)
        if ( strcmp(argv[1], ZoneID[indZ])==0 )
            break;
    if ( indZ >= NO_ZONES_ID ) {
        printf ("INV.Struct_ID = '%s'\n", argv[1]);
        return -1;
    }
    
#ifdef DEBUG_TRACE_TSLD
    sprintf(buffr, "%smutrace.txt", FOLDER_OUT_DATA);
    tsld = fopen(buffr, "w");
#endif
    
	if ( LoadDNKset(FOLDER_HUMAN, XRO_SET_FILENAME) <= 0 )  //made 'markTCX_TAG()'-->SET_TCx_TAG(Xtag)
		return -1;
    
    for ( int Z=0; Z<NO_ZONES_ID; Z++)  {
        if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[Z], 1 ) < 0 )      // made 'SET_MUZTWO_TAG'
            return -1;
    }

    if ( LoadXroZones( FOLDER_IN_DATA, ZoneID[indZ] ) < 0 )      // made 'SET_MUZONE_TAG' & not'SET_CORN_TAG'
        return -1;
    
    printf ("  calcCommnCounts()\n");
    for (int Xr=0; Xr<NO_HUMAN_XRO; Xr++ )
        vecDNK[Xr].calcCommnCounts();   //clcUnZonesArea( );
 
    sprintf(buffr, "%s%s_tst.txt", FOLDER_OUT_DATA, argv[1]);        //"%s%s_mu.txt"    for testing
    fRezult = fopen(buffr, "w");
    clock_t start;      //= clock();
    clock_t finish;     //= clock();
    double duration;    //= (double)(finish - start) / CLOCKS_PER_SEC;

    printTestCap( fRezult );                                                //    for testing
    for ( int indCancer=0; indCancer<NO_CANCER_ID; indCancer++ ) {
        if ( ! vecSAMPL.empty() )
            vecSAMPL.clear();
        if ( loadMutation( CancerID[indCancer] ) < 0 )
            return -1;
        for ( int nS=0; nS<vecSAMPL.size(); nS++ ) { 
            for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
                vecDNK[nX].calcStatMut ( vecSAMPL[nS].vMutNuc[nX] );
            }
            vecSAMPL[nS].printTestLine( fRezult, indCancer, indZ );       
        }

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
////////////////////////////////////////////////////////

void printTableCap( FILE *fRez )
{
    
//  fprintf(fRez, "Struc\tCancer\t#S\tsample\tsumIN\tsumOUTcur\tsumOUTall\tsizeIN\tsizeOUTcur\tsizeOUTall");
    fprintf(fRez, "Struc\tCancer\t#S\tsample\tsumIN\t sumOUTall\tsizeIN\t sizeOUTall\tTCxIN\tTCxOUT\tAPOin\tAPOout");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO
        fprintf(fRez, "\tX.%d_IN\tX.%d_OUT\tX.%d_sIN\tX.%d_sOUT\tX.%d_TCxIN\tX.%d_TCxOUT\tX.%d_APOin\tX.%d_APOout",
            vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
            vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSamplLine( FILE *fRez, int cancNo, int zoneNo)
{
    long sumInZon=0, sumOutCurZon=0, sumOutAllZon=0;
    long sumNucInZon=0, sumNucOutCurZon=0, sumNucOutAllZon=0;
    long sumTCxInZon=0, sumTCxOutZon=0, sumAPOinZon=0, sumAPOoutZon=0;
    
    for ( int nX=0; nX<vecDNK.size(); nX++ ) {
        sumInZon  += vecDNK[nX].cntMutInZon;
        sumOutAllZon += vecDNK[nX].cntMutOutAllZon;
        sumNucInZon += vecDNK[nX].nNucInZon;
        sumNucOutAllZon += vecDNK[nX].nNucOutAllZon;
        
        sumTCxInZon += vecDNK[nX].nTCxInZon;
        sumTCxOutZon+= vecDNK[nX].nTCxOutAllZon;
        sumAPOinZon += vecDNK[nX].cntApoInZon;
        sumAPOoutZon+= vecDNK[nX].cntApoOutAllZon;
    }

    fprintf(fRez, "%s\t%s\t%d\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld",
            ZoneID[zoneNo], CancerID[cancNo],
            (int)(this-&vecSAMPL[0]),  SampName.c_str(),
            sumInZon, sumOutAllZon, sumNucInZon, sumNucOutAllZon,
            sumTCxInZon, sumTCxOutZon, sumAPOinZon, sumAPOoutZon);
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ )
        fprintf(fRez, "\t%d\t%d\t%ld\t%ld\t%d\t%d\t%d\t%d",
                vecDNK[nX].cntMutInZon, vecDNK[nX].cntMutOutAllZon,
                vecDNK[nX].nNucInZon, vecDNK[nX].nNucOutAllZon,
                vecDNK[nX].nTCxInZon, vecDNK[nX].nTCxOutAllZon,
                vecDNK[nX].cntApoInZon, vecDNK[nX].cntApoOutAllZon );
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void printTestCap( FILE *fRez )
{
    
//  fprintf(fRez, "Struc\tCancer\t#S\tsample\tsumIN\t sumOUTall\tsizeIN\t sizeOUTall\tTCxIN\tTCxOUT\tAPOin\tAPOout");
    fprintf(fRez, "Struc\tCancer\t#S\tsample\tnMU\tsumIN\tsumINotherZ\tsumOUTall");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO
        fprintf(fRez, "\tX.%d_nMU\tX.%d_IN\tX.%d_INoterZ\tX.%d_OUT",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printTestLine( FILE *fRez, int cancNo, int zoneNo)
{
    long sumMUT=0, sumInOtherZ=0;                       // for testing
    long sumInZon=0, sumOutCurZon=0, sumOutAllZon=0;
    long sumNucInZon=0, sumNucOutCurZon=0, sumNucOutAllZon=0;
    long sumTCxInZon=0, sumTCxOutZon=0, sumAPOinZon=0, sumAPOoutZon=0;
    for ( int nX=0; nX<vecDNK.size(); nX++ ) {
        sumMUT      += vMutNuc[nX].size();                   // for testing
        sumInOtherZ += vecDNK[nX].cntMuINotheZon;            // for testing
        sumInZon  += vecDNK[nX].cntMutInZon;
        //        sumOutCurZon += vecDNK[nX].cntMutOutCurZon;
        sumOutAllZon += vecDNK[nX].cntMutOutAllZon;
        sumNucInZon += vecDNK[nX].nNucInZon;
        //        sumNucOutCurZon += vecDNK[nX].nNucOutCurZon;
        sumNucOutAllZon += vecDNK[nX].nNucOutAllZon;
        
        sumTCxInZon += vecDNK[nX].nTCxInZon;
        sumTCxOutZon+= vecDNK[nX].nTCxOutAllZon;
        sumAPOinZon += vecDNK[nX].cntApoInZon;
        sumAPOoutZon+= vecDNK[nX].cntApoOutAllZon;
    }
    
    fprintf(fRez, "%s\t%s\t%d\t%s\t%ld\t%ld\t%ld\t%ld",
            ZoneID[zoneNo], CancerID[cancNo],
            (int)(this-&vecSAMPL[0]),  SampName.c_str(),
            sumMUT, sumInZon, sumInOtherZ, sumOutAllZon);
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ )
        fprintf(fRez, "\t%ld\t%d\t%ld\t%d", (long)vMutNuc[nX].size(),
                vecDNK[nX].cntMutInZon, vecDNK[nX].cntMuINotheZon, vecDNK[nX].cntMutOutAllZon);
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////
