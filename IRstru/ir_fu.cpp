//
//  ir_fu.cpp
//  irstudy
//
//  Created by Ирина Пономарева on 23/04/2021.
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

///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////

void XROSOMA::  calcIRstat( vector< MUTANT> &vMutSet )
{
    int noCross;
    vector < MUZONE >:: iterator iterZ;
    int sizCrossZ;
    int sizTst = 1024;
    char *pTstDONE = new char[sizTst];
    memset((void *)pTstDONE, 0, sizTst);
    
    for ( int p=0; p<4; p++ )
        for ( int r=0; r<4; r++ )  {
            cntLoopMU[p][r] = 0;
            for ( int a=0; a<4; a++ )
                cntOtherMU[p][r][a] = 0;
        }
//--------------
    int LastZ = (int)vecMuZon.size() - 1;
    int curZ = 0;
    while ( curZ <= LastZ ) {
        noCross = defCrossZ( curZ, &sizCrossZ );
        if ( sizCrossZ > sizTst )   {
            delete [] pTstDONE;
            sizTst = sizCrossZ;
            pTstDONE = new char[sizTst];
        }
        memset((void *)pTstDONE, 0, sizCrossZ );
        iterZ = vecMuZon.begin() + curZ;
        for ( int n=0; n<noCross; n++ )
            procCornerC( iterZ->startPos, iterZ+n, vMutSet, pTstDONE );
        
        procCrossZones(iterZ->startPos, (iterZ+noCross-1)->stopPos,
                       vMutSet, pTstDONE );
        curZ += noCross;
    }
    
    return;
}
///////////////////////////////////////////////////////////////////////

int XROSOMA:: defCrossZ( int curZ, int *sizCrossZ )
{
    int nZones=1;
    vector <MUZONE>:: iterator itZ = vecMuZon.begin() + curZ;
    
    long begPos = itZ->startPos,
         endPos = itZ->stopPos;
    
    for ( itZ=itZ+1; itZ!=vecMuZon.end(); itZ++ ) {
        if ( itZ->startPos > endPos )
            break;
        nZones++;
        endPos = itZ->stopPos;
    }
    *sizCrossZ = (int)(endPos - begPos +1);
        
    return nZones;
}
///////////////////////////////////////////////////////////////////////

int XROSOMA:: procCornerC ( long StartPos, vector < MUZONE >:: iterator iterZ,
                         vector< MUTANT> &vMutSet, char *pDoneIt )
{
    int prev;
    long indMut; //mutPos;
    
    if ( iterZ->spacer == 0 )
        return 0;
    
    long L_corn = iterZ->startPos-1 + iterZ->repeat;
    long R_corn = L_corn + iterZ->spacer - 1;           // -1 so as it last_spacer_nuc
    if ( *(Xbody + L_corn) == 'G' )
        indMut = L_corn;             // left corner of LOOP
    else
        if ( *(Xbody + R_corn) == 'C' )
             indMut = R_corn;             // right corner of LOOP
        else
            return 0;
    
    char *pXb = Xbody + indMut;
    prev =( *pXb=='C' ) ? getNucID(*(pXb-1)) : getNucID(*(pXb+1));
    if ( prev < 0 )
        return 0;
    if ( *pXb=='G' )
        prev = cmpInd[prev];
    
    int mark = (int)(indMut - (StartPos-1));          //(iterZ->startPos-1));
    if ( *(pDoneIt+mark) )
        return 0;
    
    *(pDoneIt+mark) = 1;
    if ( ! GET_MUT_TAG(Xtag+indMut) )   {
        cntLoopMU[prev][_C] += 1;              // CornerC & NO Mutation
        return 1;
    }
    
    vector < MUTANT >:: iterator itMut =
        lower_bound( vMutSet.begin( ), vMutSet.end( ), MUTANT(indMut+1), lesser_MUT);
    
    if ( itMut->nucPos != indMut+1 )  {
        printf("ERROR!!! Inv_MUT_TAG position=%ld Not found at vMutSet\n", indMut+1);
        return 0;
    }
    if ( *pXb != itMut->nucREF ) {
        printf("ERR!!! defPinMutation():: Mismatch nucREF: Pos=%ld Xbody=%c MutREF=%c\n",
               indMut+1, *pXb, itMut->nucREF);
        return 0;
    }
    
    int mutN = getNucID(itMut->nucALT);
    if ( *pXb=='G' )
        mutN = cmpInd[mutN];
    
    cntLoopMU[prev][mutN] += 1;             // CornerC & Mutation
    
    return 1;
}

///////////////////////////////////////////////////////////////////////

void XROSOMA:: procCrossZones(long StartPos, long StopPos,
                              vector< MUTANT> &vMutSet, char *pDoneIt )
{
    vector <MUTANT>:: iterator itMut;
    vector <MUZONE>:: iterator iterZ;
    PROP_NUC cooN, compN;
    char *pXb;
    
    
    for ( long nXb=StartPos-1; nXb<=StopPos-1; nXb++ ) {
        int ixD = (int)(nXb - (StartPos-1));
        if ( *(pDoneIt+ixD) )
            continue;
        
        cooN.ini();
        compN.ini();
        pXb = Xbody + nXb;
        if ( (cooN.ref =getNucID(*pXb)) < 0 )
            continue;
        compN.ref = cmpInd[cooN.ref];
        ///
        if ( (cooN.prev =getNucID(*(pXb-1))) < 0 )
            continue;
        if ( (compN.prev =getNucID(*(pXb+1))) < 0 )
            continue;
        compN.prev = cmpInd[compN.prev];
        ///
        if ( ! GET_MUT_TAG(Xtag + nXb) )   {
            cntOtherMU[cooN.prev][cooN.ref][cooN.ref] += 1;
            cntOtherMU[compN.prev][compN.ref][compN.ref] += 1;
            *(pDoneIt+ixD) = 1;
            continue;
        }
        itMut = lower_bound( vMutSet.begin( ), vMutSet.end( ), MUTANT(nXb+1), lesser_MUT);
        if ( itMut->nucPos != nXb+1 )  {
            printf("ERROR!!! Inv_MUT_TAG position=%ld Not found at vMutSet\n", nXb+1);
            continue;
        }
        if ( *pXb != itMut->nucREF ) {
            printf("ERR!!! defPinMutation():: Mismatch nucREF: Pos=%ld Xbody=%c MutREF=%c\n",
                   nXb+1, *pXb, itMut->nucREF);
            continue;
        }
        if ( (cooN.alt =getNucID(itMut->nucALT)) < 0 )
            cooN.alt = cooN.ref;
        compN.alt = cmpInd[cooN.alt];
        
        cntOtherMU[cooN.prev][cooN.ref][cooN.alt] += 1;
        cntOtherMU[compN.prev][compN.ref][compN.alt] += 1;
        *(pDoneIt+ixD) = 1;
    }
    
    return;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


