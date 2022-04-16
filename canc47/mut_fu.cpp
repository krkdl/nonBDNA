

#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;
vector < SAMPLE >  vecSAMPL;

int getMutationPart(FILE *f_Mutas);
int skipMutationPart(FILE *f_Mutas, string &sample, int XroID);
void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

///////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 )
{
    return x1.nucPos < x2.nucPos;
}
/////////////////////////////////////////////////////////////////////////////////

int loadMutation( const char *pCancName )
{
		int ch;
		int cnt1, 
            cntMut=0;
        char buff[1024];
    
    clock_t start = clock();
    
    sprintf(buff, "%ssnv_mnv_%s.tsv", FOLDER_MUT_DATA, pCancName);
    printf("LoadMutation ('%s') ...\t", pCancName);
    
	FILE *f_Mutas=fopen( buff, "rb");
    if ( f_Mutas==NULL ) {
			printf ("\nMutations_File =""%s"" not found\n", buff);
			return -1;
    }
	
    while ( (ch=fgetc( f_Mutas )) != 0x0A && ( feof( f_Mutas )) ==0 ) {}		//skip header

		
		while ( ( cnt1=getMutationPart(f_Mutas )) > 0 )
			cntMut += cnt1;
		
		if ( cnt1 < 0 )
			return cnt1;
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    //    printf(".... dT=%5.2f sec\n", duration);
    printf(" loaded [%d] dT=%5.2f sec\n", cntMut, duration);
    
    return cntMut;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

int getMutationPart(FILE *f_Mutas)
{ 
		long int f_Pos=ftell(f_Mutas);
		char cBuffer[4096];
		char chREF, chALT;
		string sample;
		int XroID;
		int NucPos;
		int nXro, nSamp;
		int retC;  

						
	if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas) == 0 )
        return 0;
	parsMutRec( cBuffer,  &XroID, &NucPos,  &chREF, &chALT, sample);

    for ( nSamp=0; nSamp<vecSAMPL.size(); nSamp++ )
        if ( sample.compare(vecSAMPL[nSamp].SampName)==0 )
            break;
    if ( nSamp >= vecSAMPL.size() )
        vecSAMPL.push_back(SAMPLE(sample));

    nXro = findXroByID( XroID );
	if ( nXro < 0 )	{
		retC = 1 + skipMutationPart(f_Mutas, sample, XroID);
        printf("\tINV.chrID '%.16s ...' >>> Ignored %d Recs\n", cBuffer, retC);
		return	retC;
	}
   
	if ( vecSAMPL[nSamp].mutSeg[nXro].noMut > 0 ) {
        printf("\nBAD MUTATION ORDER: CHRO_ID=%d  sample=%s\n",
                vecDNK[nXro].chrNum, sample.c_str ( ));
        return (-1);
	}
    
    vecDNK[nXro].testValidDNK(NucPos, chREF, chALT );
    
    if ( vecSAMPL[nSamp].vMutNuc[nXro].empty() ||
         vecSAMPL[nSamp].vMutNuc[nXro].back().nucPos <= NucPos )
            vecSAMPL[nSamp].vMutNuc[nXro].push_back( MUTANT(NucPos, chREF, chALT ));
    else {
        vector <MUTANT>::iterator Iter =
            lower_bound( vecSAMPL[nSamp].vMutNuc[nXro].begin( ),
                         vecSAMPL[nSamp].vMutNuc[nXro].end( ), MUTANT(NucPos),  lesser_MUT);
        vecSAMPL[nSamp].vMutNuc[nXro].insert( Iter, MUTANT(NucPos, chREF, chALT ) );
    }
    
	f_Pos=ftell(f_Mutas);

	while ( 1 ) {
		if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas )==0 )
			break;
		int pars_xroNum;
		string pars_samp;
		parsMutRec( cBuffer, &pars_xroNum, &NucPos,  &chREF, &chALT, pars_samp);
		if ( pars_xroNum != XroID || pars_samp.compare(sample)!= 0 )	{		
			fseek(f_Mutas, f_Pos, SEEK_SET);
			break;
		}
		vecDNK[nXro].testValidDNK(NucPos, chREF, chALT );
        if ( vecSAMPL[nSamp].vMutNuc[nXro].empty() ||
            vecSAMPL[nSamp].vMutNuc[nXro].back().nucPos <= NucPos )
            vecSAMPL[nSamp].vMutNuc[nXro].push_back( MUTANT(NucPos, chREF, chALT ));
        else {
            vector <MUTANT>::iterator Iter =
            lower_bound( vecSAMPL[nSamp].vMutNuc[nXro].begin( ),
                        vecSAMPL[nSamp].vMutNuc[nXro].end( ), MUTANT(NucPos),  lesser_MUT);
            vecSAMPL[nSamp].vMutNuc[nXro].insert( Iter, MUTANT(NucPos, chREF, chALT ) );
        }


		f_Pos = ftell(f_Mutas);
	}
	
	return (int)vecSAMPL[nSamp].vMutNuc[nXro].size();
}
//////////////////////////////////////////////////////////

void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT, string &samp )
{

	if ( *pB == 'X' || *pB=='x' )	
		*xrNum = 23;
	else
	if ( *pB == 'Y' || *pB=='y' )
		*xrNum = 24;
	else
		*xrNum = atoi( pB );

	while ( *pB && *pB!='\t' ) pB++;	//skip CHROM
	if (*pB)	pB++;

	*Xpos = atoi( pB );								//POS
	while ( *pB && *pB!='\t' ) pB++;	//skip POS
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip ID
	if (*pB)	pB++;

	*chREF = *pB;											//REF
	while ( *pB && *pB!='\t' ) pB++;	//skip REF
	if (*pB)	pB++;

	*chALT = *pB;											//ALT
	while ( *pB && *pB!='\t' ) pB++;	//skip ALT
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip QUAL
	if (*pB)	pB++;

	while ( *pB && *pB!='\t' ) pB++;	//skip FILTER
	if (*pB)	pB++;
	
	while ( *pB && *pB!='\t' ) pB++;	//skip INFO
	if (*pB)	pB++;

	char *pBend = pB;
	while ( *pBend && *pBend!='\t' ) pBend ++;
	*pBend = '\0';
	samp = pB;

	return;
}
///////////////////////////////////////////////////////

int skipMutationPart(FILE *f_Mutas, string &sample, int XroID)
{
		long int f_Pos=ftell(f_Mutas);
		char cBuffer[4096];
		int retC = 0;

	while ( 1 ) {
        if ( fgets_ShortRec(cBuffer, sizeof(cBuffer), f_Mutas )==0 )
            break;
//		if ( fgets(cBuffer, sizeof(cBuffer), f_Mutas )==NULL )
//			break;
		int pars_xroNum;
		string pars_samp;
		int NucPos;
		char chREF, chALT;
		parsMutRec( cBuffer, &pars_xroNum, &NucPos,  &chREF, &chALT, pars_samp);
		if ( pars_xroNum != XroID || pars_samp.compare(sample)!= 0 )	{
			fseek(f_Mutas, f_Pos, SEEK_SET);
			return (retC);
		}
		retC++;
		f_Pos = ftell(f_Mutas);
	}

	return retC;
}
///////////////////////////////////////////////////////////////////////////////

void XROSOMA:: calcCommnCounts( )
{
    char *pXt;
    char *pXend=Xtag+Xsize-1;
    
    nNucOutAllZon = 0;
    nNucOutCurZon = 0;
    nTCxInZon     = 0;
    nTCxOutAllZon = 0;
    for ( pXt=Xtag; pXt<=pXend; pXt++ ) {
        if ( ! GET_MUZTWO_TAG(pXt) )
            nNucOutAllZon++;
        if ( ! GET_MUZONE_TAG(pXt) )
            nNucOutCurZon++;
        
        if ( GET_TCx_TAG(pXt) ) {
            if ( ! GET_MUZTWO_TAG(pXt) )
                nTCxOutAllZon++;
            if ( GET_MUZONE_TAG(pXt) )
                nTCxInZon++;
        }
    }

    return;
}
///////////////////////////////////////////////////////////////////////////////

void XROSOMA::  calcStatMut ( vector< MUTANT> &vMutSet )
{
    vector < MUTANT >:: iterator itMUT;
    char *pT;
    cntMutInZon     = 0;
    cntMutOutCurZon = 0;
    cntMutOutAllZon = 0;
        cntApoInZon     = 0;
        cntApoOutAllZon = 0;
    cntMuINotheZon = 0;   // for testing
    
    for ( itMUT=vMutSet.begin(); itMUT!=vMutSet.end(); itMUT++ )    {
        int itsAPO = APOtest(itMUT->nucPos, itMUT->nucREF, itMUT->nucALT);
        pT = Xtag + (itMUT->nucPos - 1);
        if ( GET_MUZONE_TAG(pT) )   {
            cntMutInZon += 1;
            if ( itsAPO )
                cntApoInZon++;
        }
        else    {
            cntMutOutCurZon += 1;
            if ( ! GET_MUZTWO_TAG(pT) ) {
                cntMutOutAllZon += 1;
                if ( itsAPO )
                    cntApoOutAllZon++;
            }
            else                    // for testing
                cntMuINotheZon++;   // for testing
        }
    }
    
    return;
    
}
///////////////////////////////////////////////////////

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

