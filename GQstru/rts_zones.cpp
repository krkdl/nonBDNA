//
//  rts_zones.cpp
//  GQstudy
//
//  Created by Ирина Пономарева on 13/03/2022.
//  Copyright © 2022 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;

int readRTSzones( const char *fRTSid )
{
    FILE *fRTS=NULL;
    char fBuff[1024];
    int cntRec=0;
    int XroNo, route;
    long zStartPos, zStopPos;
    float xParm;
    char *pStart, *pStop;
        
    printf("\nread RTS_zones ...\n" );
    clock_t start = clock();
    
    if ( ( fRTS=fopen(fRTSid, "r"))==NULL )    {
        printf("ERR_OPENop File='%s'\n", fRTSid);
        return -1;
    }
    
    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff), fRTS) == NULL )
            break;
        cntRec++;
        sscanf(fBuff, "%d\t%ld\t%ld\t%f\t%d",
                    &XroNo, &zStartPos, &zStopPos, &xParm, &route);
        if ( XroNo < 0 || XroNo > 21 )  {
            printf("INV.RTS_rec: XroNo=%d\n", XroNo);
            return -1;
        }
        if ( zStartPos <=0 || zStartPos >= zStopPos )  {
            printf("INV.RTS_rec: START:STOP = %ld : %ld\n", zStartPos, zStopPos);
            return -1;
        }

        pStart = vecDNK[XroNo].Xtag + zStartPos-1;
        pStop  = vecDNK[XroNo].Xtag + zStopPos-1;
        for ( char *pX=pStart; pX<=pStop;  pX++ )   {
            if ( route==0 )
                SET_RTS_BK0(pX);
            else
            if ( route==1 )
                SET_RTS_FW1(pX);
        }
        
    }
    printf("\treading %d zones\n", cntRec);
    return cntRec;
}
