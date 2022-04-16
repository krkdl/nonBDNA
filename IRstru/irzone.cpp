//
//  irzone.cpp
//  irstudy
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
            for ( int a=0; a<4; a++ )   {
                cntOtherMU[p][r][a][0] = 0;
                cntOtherMU[p][r][a][1] = 0;
            }
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
        procCrossZones(iterZ->startPos, (iterZ+noCross-1)->stopPos,
                       vMutSet, pTstDONE );
        curZ += noCross;
    }
    delete [] pTstDONE;
    
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

int XROSOMA:: procCornerC ( long StartPos, long indCorn, //vector < MUZONE >:: iterator iterZ,
                           vector< MUTANT> &vMutSet, char *pDoneIt )
{
    int prev;
    int mark;

    if ( ! GET_CORN_TAG(Xtag+indCorn) )
        return 0;
    
    mark = (int)(indCorn - (StartPos-1));          //(iterZ->startPos-1));
    if ( *(pDoneIt+mark) )
        return 0;
    
    char *pXb = Xbody + indCorn;
    prev =( *pXb=='C' ) ? getNucIDn(*(pXb-1)) : getNucIDn(*(pXb+1));
    if ( prev < 0 )
        return 0;
    if ( *pXb=='G' )
        prev = cmpInd[prev];
    
    *(pDoneIt+mark) = 1;
    
    if ( ! GET_MUT_TAG(Xtag+indCorn) )   {
        cntLoopMU[prev][_C] += 1;              // CornerC & NO Mutation
        return 1;
    }
    
    vector < MUTANT >:: iterator itMut =
    lower_bound( vMutSet.begin( ), vMutSet.end( ), MUTANT(indCorn+1), lesser_MUT);
    
    if ( itMut->nucPos != indCorn+1 )  {
        printf("ERROR!!! Inv_MUT_TAG position=%ld Not found at vMutSet\n", indCorn+1);
        return 0;
    }
    if ( *pXb != itMut->nucREF ) {
        printf("ERR!!! defPinMutation():: Mismatch nucREF: Pos=%ld Xbody=%c MutREF=%c\n",
               indCorn+1, *pXb, itMut->nucREF);
        return 0;
    }
    
    int mutN = getNucIDn(itMut->nucALT);
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
        if ( GET_CORN_TAG(Xtag+nXb) ) {
            procCornerC ( StartPos, (nXb), vMutSet, pDoneIt);
            continue;
        }
        cooN.ini();
        compN.ini();
        pXb = Xbody + nXb;
        if ( (cooN.ref =getNucIDn(*pXb)) < 0 )
            continue;
        compN.ref = cmpInd[cooN.ref];
        ///
        if ( (cooN.prev =getNucIDn(*(pXb-1))) < 0 )
            continue;
        if ( (compN.prev =getNucIDn(*(pXb+1))) < 0 )
            continue;
        compN.prev = cmpInd[compN.prev];
        ///
        int inLoop = GET_LOOP_TAG(Xtag + nXb);
        if ( ! GET_MUT_TAG(Xtag + nXb) )   {
            cntOtherMU[cooN.prev][cooN.ref][cooN.ref][inLoop] += 1;
            cntOtherMU[compN.prev][compN.ref][compN.ref][inLoop] += 1;
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
        
        if ( (cooN.alt =getNucIDn(itMut->nucALT)) < 0 )
            cooN.alt = cooN.ref;
        compN.alt = cmpInd[cooN.alt];
        
        cntOtherMU[cooN.prev][cooN.ref][cooN.alt][inLoop] += 1;
        cntOtherMU[compN.prev][compN.ref][compN.alt][inLoop] += 1;
        *(pDoneIt+ixD) = 1;
    }
    
    return;
}
///////////////////////////////////////////////////////////////////////

void printIRCap( FILE *fLoopRez, FILE *fRez )
{
    fprintf(fLoopRez, "#\tsample\tpC\tmut\tCnt");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
        fprintf(fLoopRez, "\tX.%d_pC\tX.%d_mut\tX.%d_Cnt",
                vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    fprintf(fLoopRez, "\n");
    
    fprintf(fRez, "#\tsample\tpC\tnuc\tmut\tLOOP\tCnt");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
        fprintf(fRez, "\tX.%d_pC\tX.%d_nuc\tX.%d_mut\tX.%d_Loop\tX.%d_Cnt",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    fprintf(fRez, "\n");
    
    return;
}
///////////////////////////////////////////////////////////////////////

void printIRsampl( FILE *fLoopRez, FILE *fRez, int nSamp )
{
    int SUMcntLoopMU[4][4];      // [prio] [alt]=[_ G A T ]
    int SUMcntOtherMU[4][4][4][2];  // [prio] [ref] [alt]=[C G A T]  [inLOOP]=[1 0]
    int prioN, refN, altN, inLoop;
    vector <XROSOMA>:: iterator itXro;
    
    for ( prioN=0; prioN<4; prioN++ )
        for ( refN=0; refN<4; refN++ )  {
            SUMcntLoopMU[prioN][refN] = 0;
            for ( altN=0; altN<4; altN++ )  {
                SUMcntOtherMU[prioN][refN][altN][0] = 0;
                SUMcntOtherMU[prioN][refN][altN][1] = 0;
            }
        }
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)  {
        itXro = vecDNK.begin()+nX;
        for ( prioN=0; prioN<4; prioN++ )
            for ( refN=0; refN<4; refN++ )  {
                SUMcntLoopMU[prioN][refN] += itXro->cntLoopMU[prioN][refN];
                for ( altN=0; altN<4; altN++ )  {
                    SUMcntOtherMU[prioN][refN][altN][0] += itXro->cntOtherMU[prioN][refN][altN][0];
                    SUMcntOtherMU[prioN][refN][altN][1] += itXro->cntOtherMU[prioN][refN][altN][1];
                }
            }
    }
    
    for ( prioN=0; prioN<4; prioN++ )
        for ( altN=0; altN<4; altN++ )  {
            fprintf(fLoopRez, "%d\t%s\t%c\t%c\t%d", nSamp, vecSAMPL[nSamp].SampName.c_str(),
                    NucID[prioN], ((altN==0) ? '-' : NucID[altN]), SUMcntLoopMU[prioN][altN]);
            for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
                fprintf(fLoopRez, "\t%c\t%c\t%d", NucID[prioN],
                        ((altN==0) ? '-' : NucID[altN]), vecDNK[nX].cntLoopMU[prioN][altN]);
            
            fprintf(fLoopRez, "\n");
        }
    
    for ( prioN=0; prioN<4; prioN++ )
        for ( refN=0; refN<4; refN++ )
            for ( altN=0; altN<4; altN++ )
                for ( inLoop=0; inLoop<2; inLoop++ )    {
                  fprintf(fRez, "%d\t%s\t%c\t%c\t%c\t%d\t%d", nSamp, vecSAMPL[nSamp].SampName.c_str(),
                        NucID[prioN], NucID[refN], ((altN==refN) ? '-' : NucID[altN]), inLoop,
                        SUMcntOtherMU[prioN][refN][altN][inLoop]);
                  for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
                    fprintf(fRez, "\t%c\t%c\t%c\t%d\t%d", NucID[prioN], NucID[refN],
                            ((altN==refN) ? '-' : NucID[altN]), inLoop,
                            vecDNK[nX].cntOtherMU[prioN][refN][altN][inLoop]);
                  fprintf(fRez, "\n");
            }
    
    return;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
