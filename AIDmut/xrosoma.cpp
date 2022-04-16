


#include <string.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
	extern FILE *tsld;
#endif

vector < XROSOMA > vecDNK;

char complment[4][2] = { {'C','G'},  {'G','C'}, {'A','T'}, {'T','A'} };
char NucID[5] = "CGAT";
int  cmpInd[4] = { _G, _C, _T, _A };

#ifdef RUN_WITH_NSBI_37
struct NSBI_37 {
    int xroID;
    long startPOS;
    long stopPOS;
} NSBI[9] = {   { 2, 75000000, 83300000}, { 3, 187439164, 187463512}, { 4, 40198526, 40246280},
                { 6, 37137921, 37143203}, { 8, 128748314, 128753679}, { 9, 36838530, 37034475},
                {10, 90740268, 90765522}, {14, 104000000, 107349540}, {22, 17900000, 25900000 },
            };
#endif

extern MUT_SIGNA   signa_C;
extern MUT_SIGNA   signa_G;
extern int curSigna;

///////////////////////////////////////////////////////

int LoadDNKset( const char *pFolder, const char *fSetName )
{
    FILE *fiSet=NULL;
    char fBuff[1024];
    int    noDNK = 0;
    long   noTcx;
    
    sprintf(fBuff, "%s%s", pFolder, fSetName);
    printf("\nLoadDNK from list '%s'.....\n", fBuff);
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
        vecDNK.push_back(XROSOMA());
        if ( (noTcx=vecDNK.back().read_DNK(pFolder, fBuff)) < 0 )    {
            vecDNK.pop_back( );
            continue;
        }
        noDNK++;
        printf("\tX.%02d [%ld]\tNoTcx=%ld\n",
               vecDNK.back().chrNum, vecDNK.back().Xsize, noTcx);
    }
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("===endLoadDNK=[%d]\tdT=%5.2f sec\n", noDNK, duration);
    
    fclose(fiSet);
    
    if ( noDNK != NO_HUMAN_XRO ) {
        printf("\n Mismatch NO_HUMAN_XRO=%d :: Must be = %d\n", noDNK, NO_HUMAN_XRO);
        return -1;
    }
    
    return noDNK;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

long XROSOMA::read_DNK( const char *pFolder, char *fName )
{
		FILE *fileDnk;
        char fPath[1024];

    sprintf(fPath, "%s%s", pFolder, fName);
	if ( !( fileDnk=fopen(fPath, "r")) )	{
		printf("File '%s': ERR_OPEN_op\n",fPath);
		return -1;
	}
	
	XfName.assign( fPath );

	fseek(fileDnk, 0, SEEK_END); 
	Xsize = ftell(fileDnk); 
    Xbody =  new char[Xsize+1];//+3];
    if ( Xsize > maxXsize )
        maxXsize = Xsize;
    

	fseek(fileDnk, 0, SEEK_SET);
		
	fgets(Xbody, Xsize, fileDnk );  // header : 

	char *pn = strstr(Xbody,"ref");
	if ( pn )
		pn = strstr(Xbody,"NC_");

	if ( !pn || Xbody[0]!='>' ) {
		printf("File '%s': not Found 1_st control Rec '%40s'\n",fPath, Xbody);
		return -1;
	}

	pn+=3;
	while( *pn && (*pn)!= '|' )  Xref_NC += *pn++;
	chrNum = atoi( Xref_NC.c_str() );	
	Xref_NC.insert(0,"NC_");

	char *pbody = Xbody;
	long count= Xsize;
	
	while ( 1 )	{
        if ( fgets(pbody, count, fileDnk) == NULL )
				break;
            if ( (pn = strchr(pbody, '\n')) )
					*pn = '\0';
			while ( *pbody ) pbody++;
			count = Xsize - (pbody-Xbody);
    }
	*(pbody+1) = '\0';
    
    Xtag  =  new char[Xsize+1];
    memset((void *)Xtag, '\0', Xsize+1);
    
#ifdef RUN_WITH_NSBI_37
    for ( int n=0; n<sizeof(NSBI); n++ )    {
        if ( NSBI[n].xroID == chrNum )  {
            StartNCBI = NSBI[n].startPOS-1; // pos from 0
            EndNCBI   = NSBI[n].stopPOS-1;  // pos from 0
            break;
        }
    }
#endif
    
    long noTcx = markTCX_TAG ( );       ///// TCx
	
    fclose(fileDnk);

	return noTcx;
}
////////////////////////////////////////////////////////////////////////////////////

int isIt_TCx( char *pXro )
{
    if ( memcmp(pXro+signa_C.offs, signa_C.sigset[curSigna], signa_C.lng)==0 )
        return 1;
    if ( memcmp(pXro+signa_G.offs, signa_G.sigset[curSigna], signa_G.lng)==0 )        
        return 2;
    
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

long XROSOMA::markTCX_TAG ( )
{
    char *pStop = Xbody + Xsize -1;
    char *pT = Xtag;
    long cnt=0;

    for ( char *pXb=Xbody; pXb<=pStop; pXb++, pT++ )  {
        if ( isIt_TCx( pXb ) ) {
            SET_TCx_TAG(pT);
            cnt++;
        }
    }
    
    return cnt;
}
////////////////////////////////////////////////////////////////////////////

int XROSOMA::testValidDNK( int Pos, char chREF, char chALT )
{

    if (  *(Xbody+Pos-1) != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%d[%d]='%c';  MUTATION: '%c' > '%c'\n",
               chrNum,Pos,*(Xbody+Pos-1), chREF, chALT);
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////

int  findXroByID( int xID )
{
    int indX;
    
    for ( indX=0; indX<vecDNK.size(); indX++ )
        if ( vecDNK[indX].chrNum == xID )
            return indX;
    printf("ERR: XroID=%d NOT FOUND\n", xID);
    return -1;
}
////////////////////////////////////////////////////// //////////////////////////////

void switchToXRO( ) //pointXRO( )
{
    for ( int n=0; n<vecDNK.size(); n++ )   {
        vecDNK[n].pTagBody = vecDNK[n].Xtag;
        vecDNK[n].TagSize  = vecDNK[n].Xsize;
    }
    return;
}
////////////////////////////////////////////////////////////////////////////////////

void switchToRT( ) //pointRT( )
{
    for ( int n=0; n<vecDNK.size(); n++ )   {
        vecDNK[n].pTagBody = vecDNK[n].RTtag;
        vecDNK[n].TagSize  = vecDNK[n].RTsize;
    }
    return;
}
////////////////////////////////////////////////////////////////////////////////////

int getNucID ( const char Nuc )
{

    switch (Nuc) {
        case 'C':
            return _C;
        case 'G':
            return _G;
        case 'A':
            return _A;
        case 'T':
            return _T;
        default:
            break;
    }
    return -1;
}
////////////////////////////////////////////////////////////////////////////////////

