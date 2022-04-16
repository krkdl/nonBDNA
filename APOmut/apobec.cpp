//
//  apobec.cpp
//  mutas
//
//  Created by Ирина Пономарева on 12/01/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;
//////////////////////////////////////////////////////////

void makeTCXbody( )
{
    char *pSelBuff = new char[XROSOMA::maxXsize+1];

    printf("\nmakeTCXbody ..... ");
    clock_t start=clock();
    
    for( int n=0; n<vecDNK.size(); n++ )    {
        memset((void *) pSelBuff, '\0', XROSOMA::maxXsize+1);
        vecDNK[n].creatTCx( pSelBuff );
    }
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\tTCXbody created : dT=%5.2f sec\n", duration);
    
    delete [] pSelBuff;
    
    return;
}
//////////////////////////////////////////////////////////

int   XROSOMA::creatTCx( char *pSelBuff )
{
    char *pSBu = pSelBuff;
    char *pXt;
    char *pXtFirst = Xtag + StartSegIndx;
    char *pXtLast  = Xtag + StopSegIndx;
    
    nTCxInZon = 0;
    nTCxOutCurZon = 0;
    nTCxOutAllZon = 0;
    for ( pXt=pXtFirst; pXt<=pXtLast; pXt++ )    {
        if ( GET_RTSEG_TAG(pXt)==0 )
            continue;
        if ( GET_TCx_TAG(pXt) ) {
            *pSBu = *pXt;
            pSBu++;
            if ( GET_MUZONE_TAG(pXt) )
                nTCxInZon++;
            else {
                nTCxOutCurZon++;
                if ( !GET_MUZTWO_TAG(pXt) )
                    nTCxOutAllZon++;
            }
        }
    }
    
    if ( TCxMemSize < (pSBu - pSelBuff) )  { //TCXsize <
        if ( TCXtag )
            delete [] TCXtag;
        TCXtag = new char [(pSBu - pSelBuff)];
        TCxMemSize = (pSBu - pSelBuff);
    }
    TCXsize = pSBu - pSelBuff;
    memcpy((void *)TCXtag, (const void *)pSelBuff, TCXsize);

    
#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\tX%d.creatTCx [%ld-%ld]=%ld\t",
            chrNum, StartSegIndx, StopSegIndx, StopSegIndx-StartSegIndx+1);
    fprintf(tsld, "RTsize=%ld inMuZo=%ld\tTCXsiz=%ld (of %ld) inMuZo=%ld\n",
            RTsize, nNucInZon, TCXsize, TCxMemSize, nTCxInZon);
#endif
    
    return (int)TCXsize;
}
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA::APOtest( long Pos, char chREF, char chALT )
{
    // APOBEC:
    //    a) TCx ---> 'C' > ['G' | 'T' ]
    //    b) xGA ---> 'G' > ['C' | 'A' ]
    // Returns: 1=APOBEC; 2=APOBEC in ZONE; 0= any other
    //
    
    char *pB = Xbody+Pos-1;
    int retC=0;
    
    if (  *pB != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum, Pos, *pB, chREF, chALT);
        return 0;
    }
    
    switch ( chREF ) {
        case 'C':
            if ( ! (chALT == 'G' || chALT == 'T') )
                return 0;
            if ( *(pB-1) != 'T' )
                return 0;
            retC = 1;
            break;
        case 'G':
            if ( ! (chALT == 'C' || chALT == 'A') )
                return 0;
            if ( *(pB+1) != 'A' )
                return 0;
            retC = 1;
            break;
        default:
            return 0;
    }
    
    return retC;
}
/////////////////////////////////////////////////////////////////

