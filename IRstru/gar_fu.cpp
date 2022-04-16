//
//  gar_fu.cpp
//  garvard
//
//  Created by Ирина Пономарева on 05/04/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"
#include "gar_fu.h"

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

extern char complment[4][2];
extern char NucID[5];
extern int  cmpInd[4];

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////

bool leseqv_PIN ( const G_PIN &x1, const G_PIN &x2 )
{
    return x1.length < x2.length;
}
///////////////////////////////////////////////////////////////////////

long findNext_C_mut (int *curMut, vector< MUTANT> &vMutSet)
{
    for ( int n=*curMut+1; n<vMutSet.size(); n++ )
        if ( vMutSet[n].nucREF == 'C' || vMutSet[n].nucREF == 'G')  {
            *curMut = n;
            return vMutSet[n].nucPos;
        }
    
    return  -1;
}
///////////////////////////////////////////////////////////////////////

int defPrioNuc(char *pXb)       //defPrioNuc(char *pXb)
{
    int nPrio;
    
    if ( ! (*pXb=='C' || *pXb=='G') )   {
        printf("ERROR!!! Nuc=%c must be 'C' or 'G'", *pXb);
        return (-1);
    }
    char Nucle = ( *pXb=='C' ) ? *(pXb-1) : *(pXb+1);
    if ( (nPrio=getNucID(Nucle)) < 0 )
        return -1;
    
    if ( *pXb=='G' )
        nPrio=cmpInd[nPrio];
    
    return nPrio;
}
///////////////////////////////////////////////////////////////////////

int defPinLen(char *pXb, char *First)
{
    char *pNuc;
    
    char *pL = (*pXb=='C') ? pXb - SPACER : pXb-1;
    char *pR = (*pXb=='C') ? pXb + 1      :  pR = pXb + SPACER;
    
    int lng =0;
    while ( pL>=First && *pR ) {
        if ( ! (pNuc=strchr(NucID, *pL)) )
            break;
        if ( complment[pNuc-NucID][1] != *pR )
            break;
        lng++;
        pL--;
        pR++;
    }
    
    return lng;
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////

void printIRCap( FILE *fLoopRez, FILE *fRez )
{
    fprintf(fLoopRez, "#\tsample\tpC\tmut\tCnt");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
        fprintf(fLoopRez, "\tX.%d_pC\tX.%d_mut\tX.%d_Cnt",
                vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    fprintf(fLoopRez, "\n");
    
    fprintf(fRez, "#\tsample\tpC\tnuc\tmut\tCnt");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
        fprintf(fRez, "\tX.%d_pC\tX.%d_nuc\tX.%d_mut\tX.%d_Cnt",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    fprintf(fRez, "\n");
    
    return;
}
///////////////////////////////////////////////////////////////////////

void printIRsampl( FILE *fLoopRez, FILE *fRez, int nSamp )
{
    int SUMcntLoopMU[4][4];      // [prio] [alt]=[_ G A T ]
    int SUMcntOtherMU[4][4][4];  // [prio] [ref] [alt]=[C G A T]
    int prioN, refN, altN;
    vector <XROSOMA>:: iterator itXro;
 
    for ( prioN=0; prioN<4; prioN++ )
        for ( refN=0; refN<4; refN++ )  {
            SUMcntLoopMU[prioN][refN] = 0;
            for ( altN=0; altN<4; altN++ )
                SUMcntOtherMU[prioN][refN][altN] = 0;
        }
 
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++)  {
        itXro = vecDNK.begin()+nX;
        for ( prioN=0; prioN<4; prioN++ )
            for ( refN=0; refN<4; refN++ )  {
                SUMcntLoopMU[prioN][refN] += itXro->cntLoopMU[prioN][refN];
                for ( altN=0; altN<4; altN++ )
                    SUMcntOtherMU[prioN][refN][altN] += itXro->cntOtherMU[prioN][refN][altN];
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
            for ( altN=0; altN<4; altN++ )  {
                fprintf(fRez, "%d\t%s\t%c\t%c\t%c\t%d", nSamp, vecSAMPL[nSamp].SampName.c_str(),
                        NucID[prioN], NucID[refN],
                        ((altN==refN) ? '-' : NucID[altN]), SUMcntOtherMU[prioN][refN][altN]);
                for ( int nX=0; nX<NO_HUMAN_XRO; nX++)
                    fprintf(fRez, "\t%c\t%c\t%c\t%d", NucID[prioN], NucID[refN],
                        ((altN==refN) ? '-' : NucID[altN]), vecDNK[nX].cntOtherMU[prioN][refN][altN]);
                fprintf(fRez, "\n");
            }

    return;
}

///////////////////////////////////////////////////////////////////////

long XROSOMA:: procCloop( vector < MUZONE >:: iterator iterZ, vector< MUTANT> &vMutSet )
{
    vector < MUTANT >:: iterator itMut;
    int prioN;
    long indMut; //mutPos;
    
    if ( iterZ->spacer == 0 )
        return -1;
    
    long L_corn = iterZ->startPos + iterZ->repeat - 1;
    long R_corn = iterZ->startPos + iterZ->repeat + iterZ->spacer - 2; // !!!!!!!!!!
    
    if ( *(Xbody + L_corn) == 'G' )
        indMut = L_corn;             // left corner of LOOP
    else
    if ( *(Xbody + R_corn) == 'C' )
        indMut = R_corn;             // right corner of LOOP
    else
        return L_corn-1;    //-1;
    
    char *pXb = Xbody + indMut;
    if ( (prioN = defPrioNuc(pXb)) < 0 )
        return -1;
    if ( ! GET_MUT_TAG(Xtag+indMut) )   {
        cntLoopMU[prioN][_C] += 1;
        return indMut;
    }
    itMut = lower_bound( vMutSet.begin( ), vMutSet.end( ), MUTANT(indMut+1), lesser_MUT);
    if ( itMut->nucPos != indMut+1 )  {
        printf("ERROR!!! Inv_MUT_TAG position=%ld Not found at vMutSet\n", indMut+1);
        return (-1);
    }
    if ( *pXb != itMut->nucREF ) {
        printf("ERR!!! defPinMutation():: Mismatch nucREF: Pos=%ld Xbody=%c MutREF=%c\n",
               indMut+1, *pXb, itMut->nucREF);
        return (-1);
    }
    int mutN = getNucID ( itMut->nucALT );
    if ( *pXb=='G' )
        mutN = cmpInd[mutN];

    cntLoopMU[prioN][mutN] += 1;
    
    return  indMut;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


