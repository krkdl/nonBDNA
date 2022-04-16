

#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

//extern FILE *fRezult;

extern vector < XROSOMA > vecDNK;
vector < SAMPLE >  vecSAMPL;

const char mutaFiles[NO_CANCER_ID][32] = {
    "snv_mnv_BLCA-US.tsv", "snv_mnv_BRCA-US.tsv", "snv_mnv_CESC-US.tsv",
    "snv_mnv_HNSC-US.tsv", "snv_mnv_LUAD-US.tsv", "snv_mnv_LUSC-US.tsv"
};

int getMutationPart(FILE *f_Mutas);//, string &samp, int *xrNum, int *nRealMutIN);
int skipMutationPart(FILE *f_Mutas, string &sample, int XroID);
void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

///////////////////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 )
{
    return x1.nucPos < x2.nucPos;
}
/////////////////////////////////////////////////////////////////////////////////

int loadMutation( int indCancer )
{
		int ch;
		int cnt1, 
            cntMut=0;
        char buff[1024];
    
    sprintf(buff, "%s%s", FOLDER_IN_DATA, mutaFiles[indCancer]);
    printf("\nLoadMutation from %s\n", buff);
    
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
			
    printf("====endLoadMutation [%d]\n", cntMut);
    
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
		int retC;   // = noMutations

						
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
		return	retC;
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
		if ( fgets(cBuffer, sizeof(cBuffer), f_Mutas )==NULL )
			break;
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
//////////////////////////////////////////////////////////////////////////////////
//              printing data in GQ structure
//              -----------------------------
//////////////////////////////////////////////////////////////////////////////////

#ifndef CALC_OUT_ZONES

void printTableCap( FILE *fRez )
{
#ifndef RTS_MODE
    fprintf(fRez, "#s\tsample\tTargStack\tTargStrnd\tAPOstack\tAPOstrnd");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) { 
        fprintf(fRez, "\tX.%d_Tstack\tX.%d_Tstrnd\tX.%d_Astack\tX.%d_Astrnd",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
#else
    fprintf(fRez, "#s\tsample\tTrgStkFw\tTrgStkBk\tTrgStk??\tTrgLinFw\tTrgLinBk\tTrgLin??");
    fprintf(fRez,           "\tApoStkFw\tApoStkBk\tApoStk??\tApoLinFw\tApoLinBk\tApoLin??");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\tX.%d_TSFw\tX.%d_TSBk\tX.%d_TS??\tX.%d_TLFw\tX.%d_TLBk\tX.%d_TL??",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
        fprintf(fRez, "\tX.%d_ASFw\tX.%d_ASBk\tX.%d_AS??\tX.%d_ALFw\tX.%d_ALBk\tX.%d_AL??",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
#endif
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSampLine( FILE *fRez)
{

    RTS_CNTS sumTarg;
        sumTarg.Stack_FW1=0; sumTarg.Stack_BK0=0; sumTarg.Stack_uDef=0;
        sumTarg.Line_FW1=0;  sumTarg.Line_BK0=0;  sumTarg.Line_uDef=0;
    RTS_CNTS sumAPO;
        sumAPO.Stack_FW1=0; sumAPO.Stack_BK0=0; sumAPO.Stack_uDef=0;
        sumAPO.Line_FW1=0;  sumAPO.Line_BK0=0;  sumAPO.Line_uDef=0;
    
    
    for ( int nX=0; nX<vecDNK.size()-2; nX++ ) { //  -2 !!!! so as noData inRTS for X_,Y_xromo
        sumTarg.Stack_FW1 += vecDNK[nX].cntTARG.Stack_FW1;
        sumTarg.Stack_BK0 += vecDNK[nX].cntTARG.Stack_BK0;
        sumTarg.Stack_uDef+= vecDNK[nX].cntTARG.Stack_uDef;
        sumTarg.Line_FW1  += vecDNK[nX].cntTARG.Line_FW1;
        sumTarg.Line_BK0  += vecDNK[nX].cntTARG.Line_BK0;
        sumTarg.Line_uDef += vecDNK[nX].cntTARG.Line_uDef;
        
        sumAPO.Stack_FW1 += cntAPO[nX].Stack_FW1;
        sumAPO.Stack_BK0 += cntAPO[nX].Stack_BK0;
        sumAPO.Stack_uDef+= cntAPO[nX].Stack_uDef;
        sumAPO.Line_FW1  += cntAPO[nX].Line_FW1;
        sumAPO.Line_BK0  += cntAPO[nX].Line_BK0;
        sumAPO.Line_uDef += cntAPO[nX].Line_uDef;
    }
    
#ifndef RTS_MODE
// fprintf(fRez, "#s\tsample\tTargStack\tTargStrnd\tAPOstack\tAPOstrnd");
    fprintf(fRez, "%d\t%s\t%d\t%d\t%d\t%d",
            (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTarg.Stack_uDef, sumTarg.Line_uDef, sumAPO.Stack_uDef, sumAPO.Line_uDef );
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%d\t%d\t%d\t%d",
                vecDNK[nX].cntTARG.Stack_uDef, vecDNK[nX].cntTARG.Line_uDef,
                cntAPO[nX].Stack_uDef, cntAPO[nX].Line_uDef);
    }
    fprintf(fRez, "\n");
#else
//  fprintf(fRez, "#s\tsample\tTrgStkFw\tTrgStkBk\tTrgStk??\tTrgLinFw\tTrgLinBk\tTrgLinB??");
//  fprintf(fRez,           "\tApoStkFw\tApoStkBk\tApoStk??\tApoLinFw\tApoLinBk\tApoLin??");
    fprintf(fRez, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
            (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTarg.Stack_FW1, sumTarg.Stack_BK0, sumTarg.Stack_uDef,
            sumTarg.Line_FW1,  sumTarg.Line_BK0,  sumTarg.Line_uDef,
            sumAPO.Stack_FW1, sumAPO.Stack_BK0, sumAPO.Stack_uDef,
            sumAPO.Line_FW1,  sumAPO.Line_BK0,  sumAPO.Line_uDef    );
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                vecDNK[nX].cntTARG.Stack_FW1, vecDNK[nX].cntTARG.Stack_BK0,
                vecDNK[nX].cntTARG.Stack_uDef,
                vecDNK[nX].cntTARG.Line_FW1, vecDNK[nX].cntTARG.Line_BK0,
                vecDNK[nX].cntTARG.Line_uDef,
                cntAPO[nX].Stack_FW1, cntAPO[nX].Stack_BK0, cntAPO[nX].Stack_uDef,
                cntAPO[nX].Line_FW1,  cntAPO[nX].Line_BK0,  cntAPO[nX].Line_uDef   );
    }
    fprintf(fRez, "\n");
#endif
    
    
    return;
}

#endif      // of ifndef CALC_OUT_ZONES
//////////////////////////////////////////////////////////////////////////////////
//              printing data for GQ structure  out ALL structures
//              --------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////

#ifdef CALC_OUT_ZONES

void printTableCap( FILE *fRez )
{
#ifndef RTS_MODE
    fprintf(fRez, "#s\tsample\tTargStack\tTargStrnd\tAPOstack\tAPOstrnd");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\tX.%d_Tstack\tX.%d_Tstrnd\tX.%d_Astack\tX.%d_Astrnd",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
#else
    fprintf(fRez, "#s\tsample\tTrgFw\tTrgBk\tTrg??");
    fprintf(fRez,           "\tApoFw\tApoBk\tApo??");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\tX.%d_T_Fw\tX.%d_T_Bk\tX.%d_T_??",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
        
        fprintf(fRez, "\tX.%d_A_Fw\tX.%d_A_Bk\tX.%d_A_??",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
#endif
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSampLine( FILE *fRez)
{
    
    RTS_CNTS sumTarg;
    sumTarg.Stack_FW1=0; sumTarg.Stack_BK0=0; sumTarg.Stack_uDef=0;
    sumTarg.Line_FW1=0;  sumTarg.Line_BK0=0;  sumTarg.Line_uDef=0;
    RTS_CNTS sumAPO;
    sumAPO.Stack_FW1=0; sumAPO.Stack_BK0=0; sumAPO.Stack_uDef=0;
    sumAPO.Line_FW1=0;  sumAPO.Line_BK0=0;  sumAPO.Line_uDef=0;
    
    
    for ( int nX=0; nX<vecDNK.size()-2; nX++ ) { //  -2 !!!! so as noData inRTS for X_,Y_xromo
        sumTarg.Stack_FW1 += vecDNK[nX].cntTARG.Stack_FW1;
        sumTarg.Stack_BK0 += vecDNK[nX].cntTARG.Stack_BK0;
        sumTarg.Stack_uDef+= vecDNK[nX].cntTARG.Stack_uDef;
        sumTarg.Line_FW1  += vecDNK[nX].cntTARG.Line_FW1;
        sumTarg.Line_BK0  += vecDNK[nX].cntTARG.Line_BK0;
        sumTarg.Line_uDef += vecDNK[nX].cntTARG.Line_uDef;
        
        sumAPO.Stack_FW1 += cntAPO[nX].Stack_FW1;
        sumAPO.Stack_BK0 += cntAPO[nX].Stack_BK0;
        sumAPO.Stack_uDef+= cntAPO[nX].Stack_uDef;
        sumAPO.Line_FW1  += cntAPO[nX].Line_FW1;
        sumAPO.Line_BK0  += cntAPO[nX].Line_BK0;
        sumAPO.Line_uDef += cntAPO[nX].Line_uDef;
    }
    
#ifndef RTS_MODE
    // fprintf(fRez, "#s\tsample\tTargStack\tTargStrnd\tAPOstack\tAPOstrnd");
    fprintf(fRez, "%d\t%s\t%d\t%d\t%d\t%d",
            (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTarg.Stack_uDef, sumTarg.Line_uDef, sumAPO.Stack_uDef, sumAPO.Line_uDef );
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%d\t%d\t%d\t%d",
                vecDNK[nX].cntTARG.Stack_uDef, vecDNK[nX].cntTARG.Line_uDef,
                cntAPO[nX].Stack_uDef, cntAPO[nX].Line_uDef);
    }
    fprintf(fRez, "\n");
#else
    //  fprintf(fRez, "#s\tsample\tTrgStkFw\tTrgStkBk\tTrgStk??\tTrgLinFw\tTrgLinBk\tTrgLinB??");
    //  fprintf(fRez,           "\tApoStkFw\tApoStkBk\tApoStk??\tApoLinFw\tApoLinBk\tApoLin??");
    fprintf(fRez, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d",
            (int)(this-&vecSAMPL[0]), SampName.c_str(),
            sumTarg.Stack_FW1, sumTarg.Stack_BK0, sumTarg.Stack_uDef,
            sumAPO.Stack_FW1, sumAPO.Stack_BK0, sumAPO.Stack_uDef
            );
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {
        fprintf(fRez, "\t%d\t%d\t%d\t%d\t%d\t%d",
                vecDNK[nX].cntTARG.Stack_FW1, vecDNK[nX].cntTARG.Stack_BK0,
                vecDNK[nX].cntTARG.Stack_uDef,
                cntAPO[nX].Stack_FW1, cntAPO[nX].Stack_BK0, cntAPO[nX].Stack_uDef
                );
    }
    fprintf(fRez, "\n");
#endif
    
    return;
}

#endif      // of defined CALC_OUT_ZONES
///////////////////////////////////////////////////////

