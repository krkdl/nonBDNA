

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
//void clearRANDStat( );
void parsMutRec( char *pB, int *xrNum, int *Xpos, char *chREF, char *chALT,
                string &samp );

//void testfseekOP (FILE *f_Mutas);

////////////////////////////////////////////////////////////////////////

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

//testfseekOP(f_Mutas);
//return -1;
		
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
	if ( vecSAMPL[nSamp].mutSeg[nXro].noMut > 0 ) {
        printf("\nBAD MUTATION ORDER: CHRO_ID=%d  sample=%s\n",
                vecDNK[nXro].chrNum, sample.c_str ( ));
        return (-1);
	}
    
    vecDNK[nXro].testValidDNK(NucPos, chREF, chALT );
    vecSAMPL[nSamp].vMutNuc[nXro].push_back( MUTANT(NucPos, chREF, chALT ));
				
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
		vecSAMPL[nSamp].vMutNuc[nXro].push_back( MUTANT(NucPos, chREF, chALT ));

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

void SAMPLE::calcREALstat( )
{
    vector < MUTANT >:: iterator itMut;
    XROSOMA *pXRO;
    int itsAPO;
    
    clearSAMPLstat();
    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
        pXRO = &vecDNK[nX];
        for ( itMut=vMutNuc[nX].begin(); itMut!=vMutNuc[nX].end(); itMut++ )    {
            char *pPos = pXRO->Xtag + itMut->nucPos-1;
            if ( ! GET_RTSEG_TAG(pPos) )
                continue;
            
            mutSeg[nX].noMut++;
            if ( (itsAPO=pXRO->APOtest(itMut->nucPos, itMut->nucREF, itMut->nucALT) ) )
                mutAPO[nX].noMut++;
                
            if ( GET_MUZONE_TAG(pPos) ) {
                mutSeg[nX].noREAL_IN++;
                if ( itsAPO )
                    mutAPO[nX].noREAL_IN++;
            }
        }
        sumMut += mutSeg[nX].noMut;
        sumMut_APO += mutAPO[nX].noMut;
        sumREAL_IN += mutSeg[nX].noREAL_IN;
        sumREAL_IN_APO += mutAPO[nX].noREAL_IN;
    }
   
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE:: correctREAL_Xmut(int nSamp)
{
    fprintf (tsld, "CorrecnMut for Sample = %d\n", nSamp);
    fprintf (tsld, "\tn_Mut\tn_IN\tn_OUTcur\tn_OUTall\tnA_Mut\tnA_IN\tnA_OUTcur\tnA_OUTall\n");
    fprintf (tsld, "sum\t%d\t%d\t \t \t%d\t%d\t\t\n",
             sumMut, sumREAL_IN,
             sumMut_APO, sumREAL_IN_APO );
    for (int nX=0; nX<NO_HUMAN_XRO; nX++ )  {
        fprintf (tsld, "X.%d\t%d\t%d\t%d\t \t%d\t%d\t%d\t \n", nX+1,
                 mutSeg[nX].noMut, mutSeg[nX].noREAL_IN, 
                 mutSeg[nX].noMut-mutSeg[nX].noREAL_IN,
                 mutAPO[nX].noMut, mutAPO[nX].noREAL_IN,
                 mutAPO[nX].noMut-mutAPO[nX].noREAL_IN);
        mutSeg[nX].noMut     -= mutAPO[nX].noMut;
        mutSeg[nX].noREAL_IN -= mutAPO[nX].noREAL_IN;
        fprintf (tsld, "--\t%d\t%d\t%d\t \n",
                mutSeg[nX].noMut, mutSeg[nX].noREAL_IN,
                mutAPO[nX].noMut-mutAPO[nX].noREAL_IN);
    }
    sumMut     -= sumMut_APO;
    sumREAL_IN -= sumREAL_IN_APO;
    fprintf (tsld, "sum--\t%d\t%d\t%d\t \n====================\n",
             sumMut, sumREAL_IN, sumMut-sumREAL_IN);
    
    return;
}
//////////////////////////////////////////////////////////////////////////////////

void SAMPLE::makeRndCycle( int nS, const char *canc_id)//f_in)
{
	int nX;

#ifdef RANDOM_TRACE
//		makeRezultName(f_in, nS, fName, "?");
//    fName += ".txt");
    sprintf(buffr, "%s%s_smp%2d.txt", FOLDER_OUT_DATA, canc_id, nS);
    FILE *fRandTrace = fopen( buffr, "w");
    
        fprintf(fRandTrace, "%s\nnoMut:\t%d\t", SampName.c_str(), sumMut_All );
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )    {
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\t%d\t", mutX[nX].noMut );
        }
        fprintf(fRandTrace, "\n");
    
        fprintf(fRandTrace, "\nseqNo\tSUM_IN\tcntLess");
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	{
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\tX.%d_INrang\t%d_cntL",
                                vecDNK[nX].chrNum, vecDNK[nX].chrNum );
        }
        fprintf(fRandTrace, "\n");
    
		fprintf(fRandTrace, "REAL\t%d\t",sumREAL_IN);
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ ) {
            if ( nX >= vecDNK.size() )
                break;
            fprintf(fRandTrace, "\t%d\t", mutX[nX].noREAL_IN);
        }
        fprintf(fRandTrace, "\n\n");
#endif
//   .........................................................
//
    printf("Samp_#%d [nM=%d - %d]...\t", nS+1, sumMut, sumMut_APO);
    clock_t start = clock();
    
    unsigned int nMu;
    unsigned long rnd_pos;
    MUTATION *pmutX, *pmutAPO;
    for ( int nCycl=0; nCycl<NO_RANDOM_CYCLES; nCycl++ )	{
        clearRANDstat( );
        for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	{
            vector <XROSOMA> :: pointer pXRO = &vecDNK[nX];
            pmutX = &mutSeg[nX];
            
            char *pPos = pXRO->pTagBody;
            long sizeBody = pXRO->TagSize;
            for ( nMu=0; nMu<pmutX->noMut; nMu++ )	{
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody;
//                if ( GET_MUZONE_TAG(pPos+rnd_pos) != 0 ) {
                if ( GET_MUZONE_TAG(pPos+rnd_pos) != 0 ) {
                    pmutX->noRAND_IN += 1;
                    sumRAND_IN += 1;
                }
            } //of noMut
            pmutAPO = &mutAPO[nX];
            pPos = pXRO->TCXtag;
            sizeBody = pXRO->TCXsize;
            for ( nMu=0; nMu<pmutAPO->noMut; nMu++ )    {
                rnd_pos = ( (double)rand() / (double)RAND_MAX ) * sizeBody; //pXRO->TCXsize;
                if ( GET_MUZONE_TAG(pPos+rnd_pos) != 0 ) {                 //pXRO->TCXtag
                    pmutAPO->noRAND_IN += 1;
                    sumRAND_IN_APO += 1;
                }
            } //of noMut_APO
            
            if ( pmutX->noREAL_IN > 0 && pmutX->noRAND_IN <= pmutX->noREAL_IN  )
                    pmutX->cntLessReal += 1;
            if ( pmutAPO->noREAL_IN > 0 && pmutAPO->noRAND_IN <= pmutAPO->noREAL_IN )
                    pmutAPO->cntLessReal += 1;
        }	// of Xromo
			
        if ( sumREAL_IN > 0 && sumRAND_IN <= sumREAL_IN )
                sumCntLess += 1;
        if ( sumREAL_IN_APO > 0 && sumRAND_IN_APO <= sumREAL_IN_APO )
                sumCntLess_APO += 1;
//========================================

#ifdef RANDOM_TRACE
			fprintf(fRandTrace, "%d\t%d\t%d", nCycl+1, sumRAND_IN, cntLessReal ); 
			for ( nX=0; nX<NO_HUMAN_XRO; nX++ )	
				if ( mutX[nX].noMut > 0 )	
					fprintf(fRandTrace, "\t%d\t%d", mutX[nX].noRAND_IN, mutX[nX].cntLessReal);
			fprintf(fRandTrace, "\n");
#endif
			
    }	// of limon cycles
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("dT=%5.2f sec\n", duration);
    
#ifdef RANDOM_TRACE
	fprintf(fRandTrace, " \t \t%.2f", (double)(cntLessReal*100) / (double)NO_RANDOM_CYCLES ); 
    for ( nX=0; nX<NO_HUMAN_XRO; nX++ )
        if ( mutX[nX].noMut > 0 )
            fprintf(fRandTrace, "\t \t%.2f", (double)(mutX[nX].cntLessReal*100) / (double)NO_RANDOM_CYCLES);
    fprintf(fRandTrace, "\n");
    
	fclose( fRandTrace );
#endif

	return;
}
////////////////////////////////////////////////////////

void printTableCap( FILE *fRez )
{
    
    fprintf(fRez, "sample\tsumSzIN\tsumMuIN\tsumSzOUT\tsumMuOUT\tSIGN");
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO
        fprintf(fRez, "\tX.%d_szIN\tX.%d_muIN\tX.%d_szOUT\tX.%d_muOUT\tX.%d_Sign",
                vecDNK[nX].chrNum, vecDNK[nX].chrNum,
                vecDNK[nX].chrNum, vecDNK[nX].chrNum, vecDNK[nX].chrNum);
    }
    fprintf(fRez, "\n");
    
    return;
}
//////////////////////////////////////////////////////////////////////////////

void SAMPLE :: printSamplLine( FILE *fRez, int APO )
{
    double dSign;
    double percIN;
    double K_in,K_out;
    long sumRTsize=0, sumRTZoneSize=0,
         sumTCXsize=0, sumTCXZoneSize=0;
    
    for ( int nX=0; nX<vecDNK.size(); nX++ ) {
        sumRTsize       += vecDNK[nX].TagSize;               // RTsize
        sumRTZoneSize   += vecDNK[nX].noNucInMuZon;          // ZONEsize in RT
        sumTCXsize      += vecDNK[nX].TCXsize;               // 
        sumTCXZoneSize  += vecDNK[nX].noTCxInMuZon; //TCXzoneSize;
    }
    
    if ( APO )
        dSign  = (double)(100 * sumCntLess_APO) / (double)NO_RANDOM_CYCLES;
    else
        dSign  = (double)(100 * sumCntLess)     / (double)NO_RANDOM_CYCLES;
//

    fprintf(fRez, "%s\t%ld\t%d\t%ld\t%d\t%.3f",
            SampName.c_str(),
            ( ( APO ) ? sumTCXZoneSize : sumRTZoneSize),
            ( ( APO ) ? sumREAL_IN_APO : sumREAL_IN ),
            ( ( APO ) ? (sumTCXsize-sumTCXZoneSize) : (sumRTsize-sumRTZoneSize)),
            ( ( APO ) ? (sumMut_APO-sumREAL_IN_APO) : (sumMut-sumREAL_IN) ),
            dSign );
    
    for ( int nX=0; nX<NO_HUMAN_XRO; nX++ ) {  //NO_HUMAN_XRO;
        if ( APO )
            dSign  = (double)(100 * mutAPO[nX].cntLessReal) / (double)NO_RANDOM_CYCLES;
        else
            dSign  = (double)(100 * mutSeg[nX].cntLessReal) / (double)NO_RANDOM_CYCLES;
        
        fprintf(fRez, "\t%ld\t%d\t%ld\t%d\t%.3f",
                ( ( APO ) ? vecDNK[nX].noTCxInMuZon : vecDNK[nX].noNucInMuZon),
                ( ( APO ) ? mutAPO[nX].noREAL_IN : mutSeg[nX].noREAL_IN),
                ( ( APO ) ? (vecDNK[nX].TCXsize-vecDNK[nX].noTCxInMuZon) : (vecDNK[nX].TagSize-vecDNK[nX].noNucInMuZon) ),
                ( ( APO ) ? (mutAPO[nX].noMut-mutAPO[nX].noREAL_IN)  : ( mutSeg[nX].noMut-mutSeg[nX].noREAL_IN) ),
                dSign);
    }
    fprintf(fRez, "\n");
    
    return;
}
////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
