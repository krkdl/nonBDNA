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

///////////////////////////////////////////////////////////////////////

int SAMPLE::  calcGQ_APO(  )
{
    vector <MUTANT>:: iterator itMU;
    char *pTg;
    char Strand;
    int Route, Atst, retC=0;

    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  { 
        cntAPO[nX].Stack_FW1 = 0; cntAPO[nX].Stack_BK0 = 0; cntAPO[nX].Stack_uDef= 0;
        cntAPO[nX].Line_FW1  = 0; cntAPO[nX].Line_BK0  = 0; cntAPO[nX].Line_uDef = 0;
        itMU = vMutNuc[nX].begin();
        for ( itMU=vMutNuc[nX].begin(); itMU !=vMutNuc[nX].end(); itMU++ ) {
            pTg = vecDNK[nX].Xtag +itMU->nucPos-1;
            if ( !( GET_TCx_TAG(pTg) ) )
                continue;
            
#ifdef CALC_OUT_ZONES
            if ( GET_MUZTWO_TAG(pTg) )      // it's in Zones
                continue;
            Strand = ( itMU->nucREF=='C' ) ? '+' : '-';
#else
            if ( ! GET_MUZTWO_TAG(pTg) )    // it's out Zones
                continue;
            Strand = ( GET_MUZONE_TAG(pTg) ) ? '+' : '-'; // (Atst>0) ? '+' : '-';
#endif
            if ( (Atst=vecDNK[nX].GQ_APOtest(itMU->nucPos, itMU->nucREF, itMU->nucALT) )==0 )
                continue;
            
            Route = -1;
            if ( GET_RTS_FW1(pTg) )
                Route = 1;
            if ( GET_RTS_BK0(pTg) )
                Route = 0;
            retC += calcRoutStrand( itMU->nucREF, Strand, Route, cntAPO[nX]);

        }
    }
    
    return retC;
}
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: GQ_APOtest( long Pos, char chREF, char chALT )
{
    
    char *pBd = Xbody+Pos-1;
//    char *pTg = Xtag +Pos-1;
    vector <MUZONE>:: iterator iterZ;
//   MUZONE mz;
//  mz.startPos = Pos;
    
    if (  *pBd != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum, Pos, *pBd, chREF, chALT);
        return 0;
    }
    
    switch ( chREF ) {
        case 'C':
            if ( *(pBd-1) != 'T' )
                return 0;
            if ( ! (chALT == 'G' || chALT == 'T') )
                return 0;

            return 1;
            
        case 'G':
            if ( *(pBd+1) != 'A' )
                return 0;
            if ( ! (chALT == 'C' || chALT == 'A') )
                return 0;

            return 1;
            
        default:
            break;
    }
       
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

int  calcRoutStrand( char Nuc, char Strand, int Route, RTS_CNTS &cnts )
{
    //                                 Zstrand  Route  RetC
    //    a) TCx :: 'C'-> ['G'|'T']     '+'     FW1     1   Stack_Forward
    //                                  '+'     BK0     2   Stack_Backward
    //                                  '+'   unKnown   5   Stack_Un
    //                                  '-'     FW1     3   Line_Forward
    //                                  '-'     BK0     4   Line_Backward
    //                                  '-'   unKnown   6   Line_Un
    //    b) xGA :: 'G'-> ['C'|'A']     '+'     FW1     4   Line_Backward
    //                                  '+'     BK0     3   Line_Forward
    //                                  '+'   unKnown   6   Line_Un
    //                                  '-'     FW1     2   Stack_Backward
    //                                  '-'     BK0     1   Stack_Forward
    //                                  '-'   unKnown   5   Stack_Un
    //
    
    switch ( Nuc ) {
        case 'C':
            switch ( Strand ) {
                case '+':
                    if ( Route==1 )
                        cnts.Stack_FW1++;
                    else
                    if ( Route==0 )
                            cnts.Stack_BK0++;
                    else    cnts.Stack_uDef++;
                    return 1;       // testing
                case '-':
                    if ( Route==1 )
                        cnts.Line_FW1++;
                    else
                    if ( Route==0 )
                            cnts.Line_BK0++;
                    else    cnts.Line_uDef++;
                    return 1;       // testing;
                default:
                    return 0;
            }

        case 'G':
            switch ( Strand ) {
                case '+':
                    if ( Route==1 )
                        cnts.Line_BK0++;
                    else
                    if ( Route==0 )
                            cnts.Line_FW1++;
                    else    cnts.Line_uDef++;
                    return 1;       // testing;
                case '-':
                    if ( Route==1 )
                        cnts.Stack_BK0++;
                    else
                    if ( Route==0 )
                            cnts.Stack_FW1++;
                    else    cnts.Stack_uDef++;
                    return 1;       // testing;
                default:
                    return 0;
            }
        default:
            break;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////

int XROSOMA:: defCrossZ( int curZ, int *sizCrossZ )
{
    int nZones=1;
    vector <MUZONE>:: iterator itZ = vecMuZon.begin() + curZ;
    char tstbuf[1024];
    
    sprintf(tstbuf,"CrossZ_X.%d #Z=%d [%ld : %ld]=%d '%c'", chrNum, curZ,
            itZ->startPos, itZ->stopPos, itZ->leng, itZ->Strand);
    
    long begPos = itZ->startPos,
    endPos = itZ->stopPos;
    
    for ( itZ=itZ+1; itZ!=vecMuZon.end(); itZ++ ) {
        if ( itZ->startPos > endPos )
            break;
        if ( nZones==1 )
            printf("%s\n",tstbuf);
        printf("\t\t   #Z=%d [%ld : %ld]=%d '%c'\t", curZ+nZones,
                itZ->startPos, itZ->stopPos, itZ->leng, itZ->Strand);
        nZones++;
        endPos = itZ->stopPos;
    }
    *sizCrossZ = (int)(endPos - begPos +1);
    
    return nZones;
}
///////////////////////////////////////////////////////////////////////

int XROSOMA:: testGQ_ggg(long StartPos, long StopPos)
{
    char *pXstart = Xbody + StartPos-1;
    char *pXstop  = Xbody + StopPos -1;
    
    if ( *pXstart != 'G' || *pXstop != 'G') {
        printf("INV.GQstruc(St/Sp) at Pos=%ld -> '%.*s...'\n",
               StartPos, (int)(StopPos-StartPos+1), pXstart);
        return 0;
    }
    
    int prioG = 1;
    int cntGGG=1, cnt=0;
    char *pX=pXstart+1;
    while ( pX <= pXstop )   {
        switch ( prioG )    {
            case 1:             // before only 1 G
                if ( *pX!='G' ) {
                    cnt++;
                    prioG = 0;
                }
                else
                    cntGGG++;

            break;
                
            default:    //case 0:
                if ( *pX=='G' )
                    prioG = 1;
                else
                    cnt++;
            break;
        }
        pX++;
    }
    if ( cntGGG < cnt ) {
        printf("INV.GQstruc(low G_cnt=%d/%d) at Pos=%ld -> '%.*s...'\n",
               cntGGG, cnt, StartPos, (int)(StopPos-StartPos+1), pXstart);
        return 0;
    }
    
    return cntGGG;
}
///////////////////////////////////////////////////////////////////////
