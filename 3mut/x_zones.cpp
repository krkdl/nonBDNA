//
//  x_zones.cpp
//  mutas
//
//  Created by Ирина Пономарева on 01/12/2020.
//  Copyright © 2020 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;

///////////////////////////////////////////////////////

int LoadXroZones( const char *pFolder, const char *fSetID )
{
    FILE *fiSet=NULL;
    char fBuff[1024];
    int    noDNK = 0;
    int    XroId, nX;
    long n;
    
    sprintf(fBuff, "%s_%s.txt", pFolder, fSetID);
    printf("\nLoad X_zones = '%s' from %s ...\n", fSetID, fBuff );
    clock_t start = clock();
    
    if ( ( fiSet=fopen(fBuff, "r"))==NULL )    {
        printf("File '%s' : ERR_OPENop\n",fBuff);
        return -1;
    }
    
    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff), fiSet) == NULL )
            break;
        if ( fBuff[0]=='/' )
            continue;
        char *pn = strchr(fBuff, '\r');
        if ( pn )
            *pn = '\0';
        else    {
            pn = strchr(fBuff, '\n');
            if ( pn )    *pn = '\0';
        }
        if ( fBuff[3]=='X')         XroId = 23;
        else
            if ( fBuff[3]=='Y')     XroId = 24;
            else
                sscanf(fBuff, " %3s%d_", (pn+100), &XroId);
        
        nX = findXroByID( XroId );
        if ( nX < 0 )   {
            printf("INV. XroID in fileName = '%s'\n", fBuff);
            return -1;
        }
        if ( (n=vecDNK[nX].read_Xzones(pFolder, fBuff)) <= 0 )
            return -1;
        noDNK++;

    }
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("===LoadXroZones\tdT=%5.2f sec\n", duration);

    fclose(fiSet);
/////////////////////////////////////////////////////
#ifdef  NO_INT_ZONESIZE
    int zFreq[NO_INT_ZONESIZE]= {0,0,0,0,0,0};
    int zFup[NO_INT_ZONESIZE] = {30, 70, 120, 500, 1000, 1000000 };
    int nF, lngF,
            lngMAX = 0;
    
    sprintf(fBuff, "%slz_%s.txt", FOLDER_OUT_DATA, fSetID);
    fiSet = fopen(fBuff, "w");
        fprintf(fiSet, "xID\tnoZones\tlngMAX");
        for ( nF=0; nF<NO_INT_ZONESIZE-1; nF++)
            fprintf(fiSet, "\t<=%d", zFup[nF]);
        fprintf(fiSet, "\t > .....\n");
    for ( nX=0; nX<vecDNK.size(); nX++) {
        for ( nF=0; nF<NO_INT_ZONESIZE; nF++)     zFreq[nF] = 0;
        lngMAX = 0;
        for ( n=0; n<vecDNK[nX].vecMuZon.size(); n++)   {
            lngF = vecDNK[nX].vecMuZon[n].second - vecDNK[nX].vecMuZon[n].first + 1;
            if ( lngMAX < lngF )
                lngMAX = lngF;
                
            for ( nF=0; nF<NO_INT_ZONESIZE-1; nF++)
                if ( lngF <= zFup[nF] ) {
                    zFreq[nF] += 1;
                    break;
                }
            if ( nF == NO_INT_ZONESIZE-1)
                zFreq[NO_INT_ZONESIZE-1] += 1;
            
        }
        fprintf(fiSet, "chr%d\t%d\t%d", vecDNK[nX].chrNum,
                (int)vecDNK[nX].vecMuZon.size(), lngMAX);
        for ( nF=0; nF<NO_INT_ZONESIZE; nF++)
            fprintf(fiSet, "\t%d", zFreq[nF]);
        fprintf(fiSet, "\n");
    }
    fclose(fiSet);
    return -1;
#endif
    
    return noDNK;
}
///////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

long int XROSOMA::read_Xzones( const char *pFolder, char *fName)
{
    FILE *fiZon=NULL;
    char buff[1024],
            *pb;
    int chrNo=-1;
    long noNucInZ=0;
    char *pStart;
    char *pStop;
    
    sprintf(buff, "%s%s", pFolder, fName);
    
    if ( !( fiZon=fopen(buff, "r")) )    {
        printf("File '%s' : ERR_OPENop\n", buff);
        return -1;
    }
    fgets_ShortRec( buff, sizeof(buff)-2, fiZon );            // header :
    while ( 1 ) {
        if ( fgets_ShortRec( buff, sizeof(buff)-2, fiZon ) == 0 )
            break;
//   testing XroID
//   -------------
        pb = strstr(buff,"chr");
        if ( !pb || pb != buff) {
            printf("INV recFormat [none 'chr']:'%s'\n", buff);
            return -1;
        }
        if ( *(pb+3)=='X')      chrNo = 23;
        else
            if ( *(pb+3)=='Y')  chrNo = 24;
            else
                                chrNo = atoi( pb+3 );
        if ( chrNo != this->chrNum ) {
            printf("MISMATH CHRid in record:'%s'\n", buff);
            return -1;
        }

        unsigned long Start=-1;
        unsigned long Stop =-1;
        sscanf( buff, "%*s\t%*s\t%*s\t%ld\t%ld\t", &Start, &Stop);
        if ( Stop >= Xsize || Start > Stop )    {
            printf("INV recFormat (Start/Stop):'%s'\n", buff);
            return -1;
        }
        vecMuZon.push_back( pair <unsigned long, unsigned long>(Start, Stop) );
        noNucInZ += ( Stop-Start+1);
        pStart = Xtag + Start - 1;
        pStop  = Xtag + Stop - 1;
        for ( char *pX=pStart; pX<=pStop;  pX++ )
            SET_MUZONE_TAG(pX);
    }
//    XzoneSize = noNucInZone;      
    fclose(fiZon);
    
    return noNucInZ;
}
///////////////////////////////////////////////////////

