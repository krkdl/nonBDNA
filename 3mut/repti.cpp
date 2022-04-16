//
//  repti.cpp
//  mutas
//
//  Created by Ирина Пономарева on 27/01/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#include "repti.hpp"
#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

//FILE *f_Out;

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

struct RT_INFILES {
    char cancID[8];
    char fname[64];
} RTfiles[NO_CANCER_ID] =   {
        {"BLCA", "RT/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed"},    //BLCA
        {"BRCA", "RT/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed"},    //BRCA
        {"CESC", "RT/wgEncodeUwRepliSeqHelas3WaveSignalRep1.mybed"},  //CESC
        {"HNSC", "RT/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed"},    //HNSC
        {"LUAD", "RT/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed"},   //LUAD
        {"LUSC", "RT/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed"}    //LUSC
    } ;

vector < pair < unsigned long, unsigned long > > segmRT;
vector < REPLYTIMESET > vecRT;


/////////////////////////////////////////////////////////////////

int loadRepTimeSet( int indCancID )
{
    unsigned long vsize, segsize;
    
    printf("\nLoadRepTimeSet for '%s'.....\n", RTfiles[indCancID].cancID);
    clock_t start = clock();

    if ( ( vsize=readRTset(indCancID) ) == 0 )
        return -1;

    segsize = vsize / NO_RT_SEGMENTS;
    printf ( "noRT=%ld SegSize=%ld\n",vsize, segsize);
    unsigned long fst=0;
    for ( int n=0; n<NO_RT_SEGMENTS-1; n++ ) {
        segmRT.push_back( pair <unsigned long, unsigned long>(fst, fst+segsize-1) );
        fst += segsize;
    }
    segmRT.push_back( pair <unsigned long, unsigned long>(fst, vsize-1) );  // It's tail
    
    sort( vecRT.begin(), vecRT.end());

    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("===end LoadRepTimeSet [%ld]\tdT=%5.2f sec\n", vsize, duration);
    
#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\n===LoadRepTimeSet for '%s'\n", RTfiles[indCancID].cancID);
    for ( int n=0; n<NO_RT_SEGMENTS; n++ )
        fprintf(tsld, "Seg.%d=[%ld-%ld]\tFirst: X%02d pos=%ld RT=%.3f\t Last: X%02d pos=%ld RT=%.3f\n",
            n, segmRT[n].first, segmRT[n].second,
            vecDNK[vecRT[segmRT[n].first].indXro].chrNum,
                vecRT[segmRT[n].first].startPos, vecRT[segmRT[n].first].repT,
            vecDNK[vecRT[segmRT[n].second].indXro].chrNum,
                vecRT[segmRT[n].second].startPos, vecRT[segmRT[n].second].repT
                );
#endif
    
    return 1;
}
/////////////////////////////////////////////////////////////////

unsigned long readRTset( int indCancID )
{
    unsigned long cntRec=0;
    FILE * RT_file;
    char buff[1024];
    int chID;
    char chr[16];
    long int Start, Stop;
    float rt;
    
    sprintf(buff, "%s%s", FOLDER_IN_DATA, RTfiles[indCancID].fname );
    RT_file = fopen(buff, "r");
    if ( ! RT_file ) {
        printf ("NOT FOUND RT_file = '%s'\n", buff);
        return 0;
    }
    
    while ( fgets(buff, sizeof(buff)-1, RT_file) )  {
        cntRec++;
        int ns = sscanf(buff, "%6s\t%ld\t%ld\t%f",
                        chr, &Start, &Stop, &rt);
        if ( ns != 4 )
            continue;
        if ( chr[3]=='X' || chr[3]=='x' )
            chID = 23;
        else
            if ( chr[3]=='Y' || chr[3]=='y' )
                chID = 24;
            else
                chID = atoi(&(chr[3]));
    
        if ( vecDNK[chID-1].chrNum == chID )
            ns = chID-1;
        else
            if ( (ns = findXroByID(chID)) < 0 ) {
                printf ("INV XroID = '%s' in rec='%s'\n", chr, buff);
                continue;
            }
        
        vecRT.push_back(REPLYTIMESET(ns,Start,Stop,rt));
        
    }
    fclose ( RT_file );
    if ( cntRec != vecRT.size() )  {
        printf ("Mismatch: read=%ld  scan=%ld\n", cntRec,vecRT.size());
        return 0;
    }
    
    return cntRec;
}
/////////////////////////////////////////////////////////////////

void markAllRTtag( )
//    mark RTSEG_TAG in Xtag  for ALL Nuc.
//    Using ay case run without RTsegmenys
{
    char *pX, *pXlast;
    
    for( int nX=0; nX<vecDNK.size(); nX++ )    {
        vecDNK[nX].StartSegIndx = 0;
        vecDNK[nX].StopSegIndx  = vecDNK[nX].Xsize-1;
        pXlast = vecDNK[nX].Xtag + vecDNK[nX].Xsize-1;
        vecDNK[nX].noNucInMuZon = 0;                  //!!!!!!
        for( pX=vecDNK[nX].Xtag; pX<=pXlast; pX++ ) {
            SET_RTSEG_TAG(pX);
            if ( GET_MUZONE_TAG(pX) )
                vecDNK[nX].noNucInMuZon++;
        }
    }
    return;
}
/////////////////////////////////////////////////////////////////

void XROSOMA:: markRTtag ( int setRT, long First, long Last )
{
    char *pStart= Xtag + First-1,
         *pStop = Xtag + Last -1,
         *pPos;
    char mask = ~(RTSEG_TAG);
    
#if defined (DEBUG_TRACE_RTSEG)
    if ( !setRT )   {
        if ( (First-StartSegIndx) < 3 )
            fprintf( tsld, "\tX%dclear [%ld", chrNum, First-1);
        if ( (StopSegIndx-Last) < 3 )
            fprintf( tsld, "\tX%dclear %ld]\t", chrNum, Last-1);
    }
#endif
    for ( pPos=pStart; pPos<=pStop; pPos++ )
        if ( setRT )
            SET_RTSEG_TAG(pPos);
        else
            *pPos &= mask;
    return;
}
/////////////////////////////////////////////////////////////////

void XROSOMA:: clearRTSEG_TAG ( int noSeg )
{
    
    if ( noSeg < 0 )
        return;
    
    long segStart=segmRT[noSeg].first,
         segStop=segmRT[noSeg].second;
    long noPos = 0;
    
#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\t\tX%d.clearRTSEG_TAG(%d)\t", chrNum, noSeg);
#endif
    
    for ( long nRT=segStart; nRT<=segStop; nRT++ )  {
        if ( vecDNK[vecRT[nRT].indXro].chrNum != chrNum )
            continue;
        markRTtag ( 0, vecRT[nRT].startPos, vecRT[nRT].stopPos );
        noPos += (vecRT[nRT].stopPos - vecRT[nRT].startPos +1);
    }
#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\t--- cleared %ld pos\n", noPos);
#endif

    return;
}
/////////////////////////////////////////////////////////////////

void makeRTSEGment( int noSeg )
//    mark RTSEG_TAG in Xtag  & create RTtag  at all XROMO
{
    char *pSelBuff = new char[XROSOMA::maxXsize+1];
    long segStart=segmRT[noSeg].first,
         segStop=segmRT[noSeg].second;
    long firstPos,
         lastPos;
    XROSOMA * pXROM;

#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\n===makeRTSEGment(%d)[%ld-%ld] RT=%.3f : %.3f\n",
            noSeg, segStart, segStop, vecRT[segStart].repT, vecRT[segStop].repT);
#endif
    printf("\nmakeRTSEGment(%d) .....\n", noSeg);
    clock_t start=clock();
    
    for( int nX=0; nX<vecDNK.size(); nX++ )    {
        memset((void *) pSelBuff, '\0', XROSOMA::maxXsize+1);
        pXROM = &vecDNK[nX];
        if ( noSeg > 0 )
            pXROM->clearRTSEG_TAG ( noSeg-1  );
        firstPos=pXROM->Xsize;
        lastPos =0;
        for ( long nRT=segStart; nRT<=segStop; nRT++ )  {
            if ( vecRT[nRT].indXro != nX )
                continue;
            pXROM->markRTtag ( 1, vecRT[nRT].startPos, vecRT[nRT].stopPos );
            if ( firstPos > vecRT[nRT].startPos )
                firstPos = vecRT[nRT].startPos;
            if ( lastPos < vecRT[nRT].stopPos )
                lastPos = vecRT[nRT].stopPos;
        }
        pXROM->creatRTSEGbody( firstPos, lastPos, pSelBuff );
    }
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("RTSEGbody created : dT=%5.2f sec\n", duration);

    delete [] pSelBuff;
    
    return;
}
//////////////////////////////////////////////////////////

int   XROSOMA::creatRTSEGbody( long firstPos, long lastPos, char *pSelBuff )
//    build tag_set 'RTSEGtag'
{
    char *pSBu = pSelBuff;
    char *pXt;
    
    if ( firstPos >= lastPos)
        return 0;
    
    StartSegIndx = firstPos - 1;
    StopSegIndx  = lastPos  - 1;
    noNucInMuZon = 0;
    char *pXtLast = Xtag + StopSegIndx;
    
    for ( pXt=Xtag+StartSegIndx; pXt<=pXtLast; pXt++ )    {
        if ( GET_RTSEG_TAG(pXt)==0 )
            continue;
        if ( GET_MUZONE_TAG(pXt) )
            noNucInMuZon++;
        *pSBu = *pXt;
        pSBu++;
    }

    if ( RTMemSize < (pSBu - pSelBuff) )  {
        if ( RTtag )
            delete [] RTtag;
        RTtag = new char [(pSBu - pSelBuff)];
        RTMemSize = (pSBu - pSelBuff);
    }
    RTsize = pSBu - pSelBuff;
    memcpy((void *)RTtag, (const void *)pSelBuff, RTsize);
    
#if defined (DEBUG_TRACE_RTSEG)
    fprintf(tsld, "\tX%d.creatRTSEGbody [%ld-%ld]\t", chrNum, StartSegIndx, StopSegIndx);
    fprintf(tsld, "LngSegm= %ld\tRTsize= %ld (of %ld)\t inMuZon= %ld\n",
            lastPos-firstPos+1, RTsize, RTMemSize, noNucInMuZon );
#endif

    return (int)RTsize;
}
///////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////

