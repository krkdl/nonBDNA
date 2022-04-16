//
//  z_zone.cpp
//  Zstudy
//
//  Created by Ирина Пономарева on 29/04/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

extern char complment[4][2];
extern char NucID[5];
extern int  cmpInd[4];
//                        C  G  A  T  x
int Rules_C[5][5] = {   { 4, 3, 3, 5, 3 },     // C
                        { 1, 0, 0, 2, 0 },     // G
                        { 1, 0, 0, 2, 0 },     // A
                        { 7, 6, 6, 8, 6 },     // T
                        { 1, 0, 0, 2, 0 }      // x
};
//                        C  G  A  T  x
int Rules_G[5][5] = {   { 0, 3, 6, 0, 0 },     // C
                        { 1, 4, 7, 1, 1 },     // G
                        { 2, 5, 8, 2, 2 },     // A
                        { 0, 3, 6, 0, 0 },     // T
                        { 0, 3, 6, 0, 0 }      // x
};

///////////////////////////////////////////////////////////////////////

void XROSOMA::  calcZ_Targets( )
{
    cTrgIN.Cc = 0.0;  cTrgIN.cC = 0.0;  cTrgIN.CT = 0.0;  cTrgIN.TC = 0.0;
    cTrgOUT.Cc = 0.0; cTrgOUT.cC = 0.0; cTrgOUT.CT = 0.0; cTrgOUT.TC = 0.0;
    
    clock_t start=clock();

    for ( long nX=0; nX<Xsize; nX++)  {
        char *pX=Xbody+nX;
        if ( *pX=='C' || *pX=='G' )
            count_CT(nX, cTrgIN, cTrgOUT);
    }
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\tX.%d\tdT=%5.2f sec\n", chrNum, duration);
    return;
}
///////////////////////////////////////////////////////////////////////

void SAMPLE::  calcZ_Mutations( )
{
    vector <MUTANT>:: iterator itMU;
   
//    clock_t start=clock();
    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
        cMutIN[nX].Cc = 0.0;  cMutIN[nX].cC = 0.0;  cMutIN[nX].CT = 0.0;  cMutIN[nX].TC = 0.0;
        cMutOUT[nX].Cc = 0.0; cMutOUT[nX].cC = 0.0; cMutOUT[nX].CT = 0.0; cMutOUT[nX].TC = 0.0;
        
        for ( itMU=vMutNuc[nX].begin(); itMU !=vMutNuc[nX].end(); itMU++ ) {
            if ( (itMU->nucREF=='C' && itMU->nucALT=='T') ||
                 (itMU->nucREF=='G' && itMU->nucALT=='A') )
              vecDNK[nX].count_CT(itMU->nucPos-1, cMutIN[nX], cMutOUT[nX]);
        }
        
    }
//    clock_t finish = clock();
//    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
//    printf("\tdT=%5.2f sec\n", duration);
}
///////////////////////////////////////////////////////////////////////

void XROSOMA::  count_CT(long indX, CNTs &cIN, CNTs &cOUT)
{
    char *pX=Xbody+indX,
         *pTag=Xtag+indX;
    int inZone, outAll;
    int nRul, nuBef, nuAft;
        
        if ( *pX=='C' ) {
            nuBef = getNucIDn ( *(pX-1) );
            nuAft = getNucIDn ( *(pX+1) );
            nRul  = Rules_C[nuBef][nuAft];
 //           printf( "%c'%c'%c -> %d\n", *(pX-1), *pX, *(pX+1), nRul );
        }
        else
            if ( *pX=='G' ) {
                nuBef = getNucIDn ( *(pX+1) );
                nuAft = getNucIDn ( *(pX-1) );
                nRul = Rules_G[nuBef][nuAft];
//                printf( "%c'%c'%c -> %d\n", *(pX+1), *pX, *(pX-1), nRul );
            }
            else    {
                printf( "INV.enter to  count_CT (%ld='%c'", indX, *pX );
                return;
            }
    
        inZone = GET_MUZONE_TAG(pTag);
        outAll = ! GET_MUZTWO_TAG(pTag);
        
        switch (nRul) {
            case 1:     // xCc  gGx
                if ( inZone )
                    cIN.Cc += 1.0;
                else
                    if ( outAll )
                        cOUT.Cc += 1.0;
                break;
                
            case 2:     // xCt  aGx
                if ( inZone )
                    cIN.CT += 1.0;
                else
                    if ( outAll )
                        cOUT.CT += 1.0;
                break;
                
            case 3:     // cCx  xGg
                if ( inZone )
                    cIN.cC += 1.0;
                else
                    if ( outAll )
                        cOUT.cC += 1.0;
                break;
                
            case 4:     // cCc  gGg
                if ( inZone )   {
                    cIN.cC += 0.5;
                    cIN.Cc += 0.5;
                }
                else
                    if ( outAll )   {
                        cOUT.cC += 0.5;
                        cOUT.Cc += 0.5;
                    }
                break;
                
            case 5:     // cCt  aGg
                if ( inZone )   {
                    cIN.cC += 0.5;
                    cIN.CT += 0.5;
                }
                else
                    if ( outAll )   {
                        cOUT.cC += 0.5;
                        cOUT.CT += 0.5;
                    }
                break;
                
            case 6:     // tCx  xGa
                if ( inZone )
                    cIN.TC += 1.0;
                else
                    if ( outAll )
                        cOUT.TC += 1.0;
                break;
                
            case 7:     // tCc  gGa
                if ( inZone )   {
                    cIN.TC += 0.5;
                    cIN.Cc += 0.5;
                }
                else
                    if ( outAll )   {
                        cOUT.TC += 0.5;
                        cOUT.Cc += 0.5;
                    }
                break;
                
            case 8:     // xCt  aGx
                if ( inZone )   {
                    cIN.TC += 0.5;
                    cIN.CT += 0.5;
                }
                else
                    if ( outAll )   {
                        cOUT.TC += 0.5;
                        cOUT.CT += 0.5;
                    }
                
                break;
            default:    // 0  = xCx, xGx
                break;
        }
    return;
}
///////////////////////////////////////////////////////////////////////
