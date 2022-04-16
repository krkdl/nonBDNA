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

extern char complment[4][2];
extern char NucID[5];
extern int  cmpInd[4];

///////////////////////////////////////////////////////

bool lesser_ZON ( const MUZONE &x1, const MUZONE &x2 )
{
    return x1.startPos < x2.startPos;
}
///////////////////////////////////////////////////////

int LoadXroZones( const char *pFolder, const char *fSetID, int TagOnly )
{
    FILE *fiSet=NULL;
    char fBuff[1024];
    int Z_size24=0;
    int    XroId, nX;
//    long n;
    
    sprintf(fBuff, "%s_%s.txt", pFolder, fSetID);
    printf("Load X_Zones('%s', %d) ...\t", fSetID, TagOnly);
    clock_t start = clock();
    
    if ( strcmp(fSetID,"cpg")==0 )  {
        Z_size24 = LoadCPGzone( FOLDER_IN_DATA, CPG_ZoneName, TagOnly);
        clock_t finish = clock();
        double duration = (double)(finish - start) / CLOCKS_PER_SEC;
        printf("Z_size=%d\tdT=%5.2f sec\n", Z_size24, duration);
        return Z_size24;
    }
    
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
        int Z_size = (int)vecDNK[nX].read_Xzones(pFolder, fBuff, TagOnly);
        if ( Z_size < 0 )
            return -1;
        if ( Z_size==0 && TagOnly==0) {
            printf("Empty ZoneArea for X.%d\n", XroId);
            return -1;
        }
        Z_size24 += Z_size;
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
    
    return Z_size24;
}
///////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

long int XROSOMA::read_Xzones( const char *pFolder, char *fName, int TagOnly )
{
    FILE *fiZon=NULL;
    char buff[1024],
            *pb;
    int chrNo=-1;
    char *pStart;
    char *pStop;
    long Z_size=0;
    int invSeq=0;
    
    sprintf(buff, "%s%s", pFolder, fName);
    
    if ( !( fiZon=fopen(buff, "r")) )    {
        printf("File '%s' : ERR_OPENop\n", buff);
        return -1;
    }
    fgets_ShortRec( buff, sizeof(buff)-2, fiZon );            // header :

    cntTARG.Stack_FW1 = 0; cntTARG.Stack_BK0 = 0; cntTARG.Stack_uDef= 0;
    cntTARG.Line_FW1  = 0; cntTARG.Line_BK0  = 0; cntTARG.Line_uDef = 0;
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

        MUZONE mz;
//                     X  Sou Type Start Stop  L Scor Strand  Rept Spac Perm
        sscanf( buff, "%*s\t%*s\t%*s\t%ld\t%ld\t%d\t%*s\t%c\t%d\t%d\t%d\t",
                        &mz.startPos, &mz.stopPos, &mz.leng, &mz.Strand, &mz.repeat, &mz.spacer, &mz.perm );
        if ( mz.stopPos >= Xsize || mz.startPos >= mz.stopPos )    {
            printf("INV recFormat (Start/Stop):'%s'\n", buff);
            return -1;
        }
        pStart = Xtag + mz.startPos - 1;
        pStop  = Xtag + mz.stopPos - 1;
        Z_size += (mz.stopPos-mz.startPos+1);
        if ( TagOnly )  {
            for ( char *pX=pStart; pX<=pStop;  pX++ )
                SET_MUZTWO_TAG(pX);
            continue;
        }
        //-----------------------------------------------------
        
        invSeq += testSeqNuc(mz, buff);
        
        vector <MUZONE>:: iterator itMUZ =vecMuZon.end();
        
        if ( vecMuZon.empty() || vecMuZon.back().startPos < mz.startPos ) {
            vecMuZon.push_back( mz );
        }
        else {
            itMUZ = lower_bound( vecMuZon.begin( ), vecMuZon.end( ), mz, lesser_ZON);
            vecMuZon.insert( itMUZ, mz );
        }

        for ( char *pX=pStart; pX<=pStop;  pX++ )   {
            if ( mz.Strand=='+' )
                SET_MUZONE_TAG(pX);
            else
                SET_MUZMNS_TAG(pX);
        }

#ifndef CALC_OUT_ZONES
        markTCX_TAG ( mz );       // move it to readDNK
#endif                                  // to make possible counting TCx out of all Zones
    }
  
    fclose(fiZon);
    if ( invSeq > 0 )   {
        printf ( "\n!!! ERROR SEQUENCE at X.%d.  See 'mutrace.txt'\n", this->chrNum);
    }
    
    return Z_size;
}
///////////////////////////////////////////////////////////////////////////////

int LoadCPGzone( const char *pFolder, const char *cpgZoneName, int TagOnly)
{
    FILE *fiZon=NULL;
    char buff[1024], *pb;
    int indX=0, chrN=-1;
    char *pStart;
    char *pStop;
    int Z_size=0;      // summ (stop-start) along all zones
    
    sprintf(buff, "%s%s", pFolder, cpgZoneName);
    if ( !( fiZon=fopen(buff, "r")) )    {
        printf("File '%s' : ERR_OPENop\n", buff);
        return -1;
    }
    
    while ( 1 ) {
        if ( fgets_ShortRec( buff, sizeof(buff)-2, fiZon ) == 0 )
            break;
        //   testing XroID
        //   -------------
        pb = strstr(buff,"chr");
        if ( !pb || pb != buff) {
            printf("INV recFormat [none 'chr']:'%s'\n", buff);
            continue;
        }
        if ( *(pb+3)=='X')      chrN = 23;
        else
            if ( *(pb+3)=='Y')  chrN = 24;
            else
                chrN = atoi( pb+3 );
        if ( chrN != vecDNK[indX].chrNum ) {
            if ( (indX = findXroByID(chrN)) < 0 )   {
                printf("\tIgnore Rec:'%.40s'...\n", buff);
                continue;
            }
        }
        
        MUZONE mz;  // set it at init ={0,0,0,0,0,0};
        //                     X  Sou Type Start Stop  L  Scor Str  Rept Spac Perm
        if ( sscanf(buff, "%*s\t%ld\t%ld\t", &mz.startPos, &mz.stopPos) != 2 ||
            mz.stopPos >= vecDNK[indX].Xsize || mz.startPos >= mz.stopPos )    {
            printf("\tINV RecFormat (Start/Stop):'%.32s'... ignored\n", buff);
            continue;
        }
        
        pStart = vecDNK[indX].Xtag + mz.startPos - 1;
        pStop  = vecDNK[indX].Xtag + mz.stopPos - 1;
        Z_size += (mz.stopPos-mz.startPos+1);
        if ( TagOnly )  {
            for ( char *pX=pStart; pX<=pStop;  pX++ )
                SET_MUZTWO_TAG(pX);
            continue;
        }
        //-----------------------------------------------------
        vector <MUZONE>:: iterator itMUZ = vecDNK[indX].vecMuZon.end();
        
        if ( vecDNK[indX].vecMuZon.empty() ||
            vecDNK[indX].vecMuZon.back().startPos < mz.startPos ) {
            vecDNK[indX].vecMuZon.push_back( mz );
        }
        else {
            itMUZ = lower_bound( vecDNK[indX].vecMuZon.begin( ),
                                vecDNK[indX].vecMuZon.end( ), mz, lesser_ZON);
            vecDNK[indX].vecMuZon.insert( itMUZ, mz );
        }
        for ( char *pX=pStart; pX<=pStop;  pX++ )   {
            if ( GET_MUZONE_TAG(pX) )       // cross_part
                continue;
            SET_MUZONE_TAG(pX);
        }
    }
    
    fclose(fiZon);
    
    return Z_size;
}
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: testSeqNuc(MUZONE &mz, char *buff)
{
    int ntab=0;
    char *pb = buff;
    
    while ( ntab < 13 ) {
        if ( *pb == '\t' )
            ntab++;
        pb++;
    }
    if ( !( mz.Strand=='+' || mz.Strand=='-') ) {
        printf("!!!Zone_X.%d Pos=%ld; INV.Strand = '%c'\n", chrNum, mz.startPos, mz.Strand);
        return 0;
    }
    
    char *pXr = Xbody + mz.startPos-1;
    int len = ( mz.leng > 700 ) ? 700 : mz.leng;
    if ( mz.Strand=='+' )   {
        for ( int n=0; n<len; n++)
            if ( toupper(*(pXr+n)) != toupper(*(pb+n)) )  {
                fprintf(tsld, "INV.ZoneSeq_X.%d Pos=%ld Xro=%.30s\n", chrNum, mz.startPos, pXr);
                return 1;
            }
    }
    if ( mz.Strand=='-' )   {
        pXr = Xbody + mz.startPos-1 + mz.leng-1;
        for ( int n=0; n<len; n++)  {
            int nid = getNucIDn( toupper(*(pb+n)) );
            if ( *(pXr-n) != NucID[cmpInd[nid]] )  {
                fprintf(tsld, "INV.ZoneSeq_X.%d Pos=%ld Xro=%.30s\n", chrNum, mz.startPos, pXr);
                return 1;
            }
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

void XROSOMA::  testZoneCross(  )
{
    int noCross=0;
    int LastZ = (int)vecMuZon.size() - 1;
    int curZ = 0;
    while ( curZ <= LastZ ) {
        noCross = getCrossZ(curZ);
        curZ += noCross;
        
    }
}
///////////////////////////////////////////////////////////////////////

int XROSOMA:: getCrossZ( int curZ )
{
    int nZones=1;
    int n;
    vector <MUZONE>:: iterator itZ = vecMuZon.begin() + curZ;
    long begPos = itZ->startPos,
        endPos = itZ->stopPos;
    
    for ( itZ=itZ+1; itZ!=vecMuZon.end(); itZ++ ) {
        if ( itZ->startPos > endPos )
            break;
        nZones++;
        endPos = itZ->stopPos;
    }
    if ( nZones==1 )
        return 1;
    
    int lng = (int)(endPos - begPos +1);
    fprintf(tsld, "X.%02d CrossZ=%d [%ld : %ld] ='%.*s'\n", chrNum, nZones,
            vecMuZon[curZ].startPos, vecMuZon[curZ+nZones-1].stopPos, lng, Xbody+begPos-1);
    
    for ( int nZ=0; nZ<nZones; nZ++ )  {
        itZ = vecMuZon.begin() + curZ + nZ;
        long np;
        fprintf(tsld, "\t'%c'#Z=%d L%d[%ld : %ld]='", itZ->Strand, curZ+nZ,
            itZ->leng, itZ->startPos, itZ->stopPos);
        for ( np=begPos-1; np<itZ->startPos-1; np++ )
            fprintf(tsld, ".");
//            fprintf(tsld, "%c", tolower(*(Xbody+np)));
        
        if ( itZ->Strand=='+' )
            for ( np=itZ->startPos-1; np<=itZ->stopPos-1; np++ )
                fprintf(tsld, "%c", *(Xbody+np));
        else
            for ( np=itZ->stopPos-1; np>=itZ->startPos-1; np-- )
                fprintf(tsld, "%c", getCompNucID (*(Xbody+np)) );
        
        for ( np=itZ->stopPos; np<endPos; np++ )
            fprintf(tsld, ".");
//            fprintf(tsld, "%c", tolower(*(Xbody+np)));
        fprintf(tsld, "'\n");
    }
        
    return nZones;
}
///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

