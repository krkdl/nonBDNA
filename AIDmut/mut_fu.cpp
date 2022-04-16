

#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif


extern vector < XROSOMA > vecDNK;
vector < SAMPLE >  vecSAMPL;

const char mutaFiles[NO_CANCER_ID][32] = {
    "snv_mnv_BLCA-US.tsv", "snv_mnv_BRCA-US.tsv", "snv_mnv_CESC-US.tsv",
    "snv_mnv_HNSC-US.tsv", "snv_mnv_LUAD-US.tsv", "snv_mnv_LUSC-US.tsv"
};

extern MUT_SIGNA   signa_C;
extern MUT_SIGNA   signa_G;
extern int curSigna;

int getMutationPart(FILE *f_Mutas);//, string &samp, int *xrNum, int *nRealMutIN);
int skipMutationPart(FILE *f_Mutas, string &sample, int XroID);
void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

////////////////////////////////////////////////////////////////////////

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
        clock_t start = clock();
    
    sprintf(buff, "%s%s", FOLDER_IN_DATA, mutaFiles[indCancer]);
    printf("\nLoadMutation from %s\n", mutaFiles[indCancer]);
    
	FILE *f_Mutas=fopen( buff, "rb");
    if ( f_Mutas==NULL ) {
			printf ("\nMutations_File =""%s"" not found\n", buff);
			return -1;
    }
	
    while ( (ch=fgetc( f_Mutas )) != 0x0A && ( feof( f_Mutas )) ==0 )
        { }	//skip header
		
		while ( ( cnt1=getMutationPart(f_Mutas )) > 0 )
			cntMut += cnt1;
		
		if ( cnt1 < 0 )
			return cnt1;
			
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\tloaded [%d] dT=%5.2f sec\n", cntMut, duration);

    return cntMut;
}
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
		printf("\tINV.chrID '%.16s ...' >>> Ignored %d Recs\n", cBuffer, retC);
        return	retC;
	}

    if ( vecSAMPL[nSamp].vMutNuc[nXro].size() > 0 ) {
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
                    vecSAMPL[nSamp].vMutNuc[nXro].end( ),
                    MUTANT(NucPos),  lesser_MUT);
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
                        vecSAMPL[nSamp].vMutNuc[nXro].end( ),
                        MUTANT(NucPos),  lesser_MUT);
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

void SAMPLE::clearSAMPLstat() {              // call first at each RT circl
    sumXr.nMut=0;         sumXrAPO.nMut=0;
    sumXr.nREAL_IN=0;     sumXrAPO.nREAL_IN=0;
    sumXr.nREAL_OUTcur=0; sumXrAPO.nREAL_OUTcur=0;
    sumXr.nREAL_OUTall=0; sumXrAPO.nREAL_OUTall=0;
    
    sumXr.nRAND_IN=0;     sumXrAPO.nRAND_IN=0;
    sumXr.nRAND_OUTcur=0; sumXrAPO.nRAND_OUTcur=0;
    sumXr.nRAND_OUTall=0; sumXrAPO.nRAND_OUTall=0;
    sumXr.cntLessReal=0;  sumXrAPO.cntLessReal=0;
    for (int n=0; n<NO_HUMAN_XRO; n++) {
        muXr[n].nMut=0;         muXrAPO[n].nMut=0;
        muXr[n].nREAL_IN=0;     muXrAPO[n].nREAL_IN=0;
        muXr[n].nREAL_OUTcur=0; muXrAPO[n].nREAL_OUTcur=0;
        muXr[n].nREAL_OUTall=0; muXrAPO[n].nREAL_OUTall=0;
        
        muXr[n].nRAND_IN=0;     muXrAPO[n].nRAND_IN=0;
        muXr[n].nRAND_OUTcur=0; muXrAPO[n].nRAND_OUTcur=0;
        muXr[n].nRAND_OUTall=0; muXrAPO[n].nRAND_OUTall=0;
        muXr[n].cntLessReal=0;  muXrAPO[n].cntLessReal=0;
    };
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE::clearRANDstat() {              // call first within each RANDOM circl
    sumXr.nRAND_IN=0;     sumXrAPO.nRAND_IN=0;
    sumXr.nRAND_OUTcur=0; sumXrAPO.nRAND_OUTcur=0;
    sumXr.nRAND_OUTall=0; sumXrAPO.nRAND_OUTall=0;
    for (int n=0; n<NO_HUMAN_XRO; n++)  {
        muXr[n].nRAND_IN=0;     muXrAPO[n].nRAND_IN=0;
        muXr[n].nRAND_OUTcur=0; muXrAPO[n].nRAND_OUTcur=0;
        muXr[n].nRAND_OUTall=0; muXrAPO[n].nRAND_OUTall=0;
    };
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE::calcREALstat( )
{
    vector < MUTANT >:: iterator itMut;
    XROSOMA *pXRO;
    int itsAPO;
    
    clearSAMPLstat();
    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
        pXRO = &vecDNK[nX];
        for ( itMut=vMutNuc[nX].begin(); itMut!=vMutNuc[nX].end(); itMut++ )    {
            long PosNo = itMut->nucPos-1;
#ifdef RUN_WITH_NSBI_37
            if ( PosNo >= pXRO->StartNCBI && PosNo <= pXRO->EndNCBI )
                continue;
#endif
            char *pPos = pXRO->Xtag + PosNo;
            
            muXr[nX].nMut++;
            if ( (itsAPO=pXRO->APOtest(itMut->nucPos, itMut->nucREF, itMut->nucALT) ) )
                muXrAPO[nX].nMut++;
                
            if ( GET_MUZONE_TAG(pPos) ) {
                muXr[nX].nREAL_IN++;
                if ( itsAPO )
                    muXrAPO[nX].nREAL_IN++;
            }
            else {
                muXr[nX].nREAL_OUTcur++;
                if ( itsAPO )
                    muXrAPO[nX].nREAL_OUTcur++;
                if ( ! GET_MUZTWO_TAG(pPos) ) {
                    muXr[nX].nREAL_OUTall++;
                    if ( itsAPO )
                        muXrAPO[nX].nREAL_OUTall++;
                }
            }
        }
        sumXr.nMut        += muXr[nX].nMut;
        sumXr.nREAL_IN    += muXr[nX].nREAL_IN;
        sumXr.nREAL_OUTcur+= muXr[nX].nREAL_OUTcur;
        sumXr.nREAL_OUTall+= muXr[nX].nREAL_OUTall;
        
        sumXrAPO.nMut        += muXrAPO[nX].nMut;
        sumXrAPO.nREAL_IN    += muXrAPO[nX].nREAL_IN;
        sumXrAPO.nREAL_OUTcur+= muXrAPO[nX].nREAL_OUTcur;
        sumXrAPO.nREAL_OUTall+= muXrAPO[nX].nREAL_OUTall;
        
    }
   
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE:: correctREAL_Xmut(int nSam)
{

    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
        muXr[nX].nMut         -= muXrAPO[nX].nMut;
        muXr[nX].nREAL_IN     -= muXrAPO[nX].nREAL_IN;
        muXr[nX].nREAL_OUTcur -= muXrAPO[nX].nREAL_OUTcur;
        muXr[nX].nREAL_OUTall -= muXrAPO[nX].nREAL_OUTall;
    }
    sumXr.nMut         -= sumXrAPO.nMut;
    sumXr.nREAL_IN     -= sumXrAPO.nREAL_IN;
    sumXr.nREAL_OUTcur -= sumXrAPO.nREAL_OUTcur;
    sumXr.nREAL_OUTall -= sumXrAPO.nREAL_OUTall;

    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE::makeRndCycle( int nS, const char *canc_id)//f_in)
{
	int nX;
    unsigned int nM;
    unsigned long rnd_pos;
    MUTATION *pmutX, *pmutAPO;
    
#if ! defined ( RUN_WITH_RANDOM_CYCLES )
    clearRANDstat( );
#endif
    
    for ( int nCycl=0; nCycl<NO_RANDOM_CYCLES; nCycl++ )	{
        clearRANDstat( );
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	{
            vector <XROSOMA> :: pointer pXRO = &vecDNK[nX];
            
            pmutX = &muXr[nX];
            char *pPos = pXRO->pTagBody;
            long sizeBody = pXRO->TagSize;
            for ( nM=0; nM<pmutX->nMut; nM++ )	{
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody;
                if ( GET_MUZONE_TAG(pPos+rnd_pos) )
                    pmutX->nRAND_IN += 1;
                else {
                    pmutX->nRAND_OUTcur++;
                    if ( ! GET_MUZTWO_TAG(pPos+rnd_pos) )
                        pmutX->nRAND_OUTall++;
                }
            } //of noMut
            sumXr.nRAND_IN     += pmutX->nRAND_IN;
            sumXr.nRAND_OUTcur += pmutX->nRAND_OUTcur;
            sumXr.nRAND_OUTall += pmutX->nRAND_OUTall;
            
            pmutAPO = &muXrAPO[nX];
            pPos = pXRO->TCXtag;
            sizeBody = pXRO->TCXsize;
            for ( nM=0; nM<pmutAPO->nMut; nM++ )    {
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody;
                if ( GET_MUZONE_TAG(pPos+rnd_pos) )
                    pmutAPO->nRAND_IN += 1;
                else {
                    pmutAPO->nRAND_OUTcur++;
                    if ( ! GET_MUZTWO_TAG(pPos+rnd_pos) )
                        pmutAPO->nRAND_OUTall++;
                }
            } //of noMut_APO
            sumXrAPO.nRAND_IN     += pmutAPO->nRAND_IN;
            sumXrAPO.nRAND_OUTcur += pmutAPO->nRAND_OUTcur;
            sumXrAPO.nRAND_OUTall += pmutAPO->nRAND_OUTall;
            
            if ( pmutX->nREAL_IN >0 && pmutX->nRAND_IN <= pmutX->nREAL_IN )
                    pmutX->cntLessReal += 1;
            if ( pmutAPO->nREAL_IN >0 && pmutAPO->nRAND_IN <= pmutAPO->nREAL_IN )
                    pmutAPO->cntLessReal += 1;
        }	// of Xromo
			
        if ( sumXr.nREAL_IN >0 && sumXr.nRAND_IN <= sumXr.nREAL_IN )
                sumXr.cntLessReal += 1;
        if ( sumXrAPO.nREAL_IN >0 && sumXrAPO.nRAND_IN <= sumXrAPO.nREAL_IN )
                sumXrAPO.cntLessReal += 1;
//========================================
	
    }	// of limon cycles
    
	return;
}
////////////////////////////////////////////////////////

void printTableCap( FILE *fRez )
{
    
    fprintf(fRez, "signat\tsample\tszIN\tmuIN\tszOUall\tmuOUall\tszOUcur\tmuOUcur\tSIGN");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO
        fprintf(fRez, "\tX.%d_szIN\tX.%d_muIN\tX.%d_szOUall\tX.%d_muOUall\tX.%d_szOUcur\tX.%d_muOUcur\tX.%d_Sign",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSamplLine( FILE *fRez, int APO )
{
    double dSign;
//    double percIN;
//    double K_in,K_out;
    long sumNucIN=0, sumTCxIN=0,
         sumNucOUtcur=0, sumTCxOUtcur=0,
         sumNucOUtall=0, sumTCxOUtall=0;
    
    for ( int nX=0; nX<vecDNK.size(); nX++ ) {
        sumNucIN   += vecDNK[nX].nNucInZon;
        sumNucOUtcur += vecDNK[nX].nNucOutCurZon;
        sumNucOUtall += vecDNK[nX].nNucOutAllZon;
        
        sumTCxIN  += vecDNK[nX].nTCxInZon;
        sumTCxOUtcur += vecDNK[nX].nTCxOutCurZon;
        sumTCxOUtall += vecDNK[nX].nTCxOutAllZon;
        
    }
    if ( NO_RANDOM_CYCLES==0)
        dSign  = 999;
    else
        dSign  = ( APO )
            ? (double)(100 * sumXrAPO.cntLessReal) / (double)NO_RANDOM_CYCLES
            : (double)(100 * sumXr.cntLessReal)    / (double)NO_RANDOM_CYCLES;
//
// fprintf(fRez, "sample\tSzIN\tMuIN\tSzOUTall\tMuOUTall\tSzOUTcur\tMuOUTcur\tSIGN");
    fprintf(fRez, "%s+%s\t%s\t%ld\t%d\t%ld\t%d\t%ld\t%d\t%.3f",
            signa_C.sigset[curSigna], signa_G.sigset[curSigna],  SampName.c_str(),
            ( ( APO ) ? sumTCxIN : sumNucIN),
            ( ( APO ) ? sumXrAPO.nREAL_IN : sumXr.nREAL_IN ),
            ( ( APO ) ? sumTCxOUtall : sumNucOUtall),
            ( ( APO ) ? sumXrAPO.nREAL_OUTall : sumXr.nREAL_OUTall ),
            ( ( APO ) ? sumTCxOUtcur : sumNucOUtcur),
            ( ( APO ) ? sumXrAPO.nREAL_OUTcur : sumXr.nREAL_OUTcur ),
            dSign );
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO;
        if ( NO_RANDOM_CYCLES==0)
            dSign  = 999;
        else
            dSign = ( APO )
                ? (double)(100 * muXrAPO[nX].cntLessReal) / (double)NO_RANDOM_CYCLES
                : (double)(100 * muXr[nX].cntLessReal)    / (double)NO_RANDOM_CYCLES;
        
//fprintf(fRez, "\tX.%d_szIN\tX.%d_muIN\tX.%d_szOUT\tX.%d_muOUT\tX.%d_Sign",
        fprintf(fRez, "\t%ld\t%d\t%ld\t%d\t%ld\t%d\t%.3f",
                ( ( APO ) ? vecDNK[nX].nTCxInZon : vecDNK[nX].nNucInZon),
                ( ( APO ) ? muXrAPO[nX].nREAL_IN : muXr[nX].nREAL_IN),
                ( ( APO ) ? vecDNK[nX].nTCxOutAllZon : vecDNK[nX].nNucOutAllZon),
                ( ( APO ) ? muXrAPO[nX].nREAL_OUTall : muXr[nX].nREAL_OUTall),
                ( ( APO ) ? vecDNK[nX].nTCxOutCurZon : vecDNK[nX].nNucOutCurZon),
                ( ( APO ) ? muXrAPO[nX].nREAL_OUTcur : muXr[nX].nREAL_OUTcur),
                dSign);
    }
    fprintf(fRez, "\n");
    
    return;
}
////////////////////////////////////////////////////////

