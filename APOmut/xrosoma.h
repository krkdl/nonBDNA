#ifndef _XROSOMA_H__
#define _XROSOMA_H__


#include "vector"
#include "string"
//#include "dnk.h"
using namespace std;

#define _C 0
#define _G 1
#define _A 2
#define _T 3

class HUGEN {
public:
    long startpos;
    long endpos;
    mutable long maxendpos;
    char strand;
    long geneID;
    string geneName;
		int count_Pos;			// ONLY for testing GetGenesByPos()

    HUGEN( long s_pos ) { startpos = s_pos; };
    HUGEN(long s_pos, long e_pos, char s_and, long g_id, string g_name)
				{ startpos=s_pos; endpos=e_pos; maxendpos=0; strand=s_and; geneID=g_id; 
					geneName.assign( g_name ); count_Pos=0; 
				};
		
};
//////////////////////////////////////////////////////////

struct MUTANT
{
		long nucPos;
		char nucREF;
		char nucALT;
		MUTANT() {nucPos=-1; nucREF='?'; nucALT='?';};
        MUTANT(long P) {nucPos=P; };
		MUTANT(long P, char R, char A) {nucPos=P; nucREF=R; nucALT=A;};
};

class MUTATION { 
public:
//    int chrNum;
			int nMut;          // Number REAL mutation in cur_RTseg
			int nREAL_IN;      //   part of 'noMut' located in cur_Zone
            int nREAL_OUTcur;  //   part of 'noMut' located outside cur_Zone
            int nREAL_OUTall;  //     "        "       "       "    all Zones
			int nRAND_IN;      // as 'nREAL_IN'     for  RANDOM case
            int nRAND_OUTcur;  // as 'nREAL_OUTcur' "       "     "
            int nRAND_OUTall;  // as 'nREAL_OUTall' "       "     "

			int cntLessReal;        //  counts where (sumRAND_IN <= sumREAL_IN)
    
    MUTATION() {  nMut=0; nREAL_IN=0; nREAL_OUTcur=0; nREAL_OUTall=0;
                          nRAND_IN=0; nRAND_OUTcur=0; nRAND_OUTall=0;
                  cntLessReal=0; };
};
bool lesser_MUT ( const MUTANT &x1, const MUTANT &x2 );

struct MUZONE {
    long startPos;
    long stopPos;
    int leng;
    int repeat;
    int spacer;
    int perm;
};
bool lesser_ZON ( const MUZONE &x1, const MUZONE &x2 );

//////////////////////////////////////////////////////////

class SAMPLE {
public:
		string SampName;

    vector < MUTANT > vMutNuc[NO_HUMAN_XRO]; // common set mutarions in XRO's
		MUTATION muXr[NO_HUMAN_XRO];  
        MUTATION sumXr;
        MUTATION muXrAPO[NO_HUMAN_XRO];  
        MUTATION sumXrAPO;
   
        SAMPLE(string S) { SampName=S;
                sumXr.nMut=0; sumXr.nREAL_IN=0; sumXr.nREAL_OUTcur=0;
                sumXr.nREAL_OUTall=0; sumXr.cntLessReal=0;
            sumXrAPO.nMut=0; sumXrAPO.nREAL_IN=0; sumXrAPO.nREAL_OUTcur=0;
            sumXrAPO.nREAL_OUTall=0; sumXrAPO.cntLessReal=0;
        };
        void clearSAMPLstat( );
		void clearRANDstat( );
        void markTCXtag(int nXro);
        void calcREALstat( );
        void correctREAL_Xmut(int nSam);
		void makeRndCycle( int nS, const char *canc_id);
        void printSamplLine( FILE *fRezint, int APO );
        void markMutations ( int ON );
};

void printTableCap( FILE *fRez );

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

class  XROSOMA {
private:
		vector <HUGEN> :: iterator begSearch;
public:
    static long maxXsize;
	string XfName;
	string Xref_NC;
	int  chrNum;				// get from Xref_NC
	string version;
//
//  ===== real body of XRO will be read from files =============================
    long Xsize;
    char *Xbody;
    char *Xtag;     // 0x01=RTsegment; 0x80=inZone; 0x40=TCx; 0x20=APOBEC;
                    // RTsegment defines XRO work_area for creating TCXtag, noMut, REAL_IN, ets.
    long StartSegIndx;      // Bounds of Current RTsegment in 'Xtag' (from 0 ....)
    long StopSegIndx;       // (from 0 ....)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  ===== Current RT segment. Will be updated for each segment. ================
    long RTMemSize;     // = sizeof (RTtag)
    char *RTtag;        // contains tags from '*Xtag' with RTSEG_TAG
    long RTsize;        // = real Segment_Size  (<=RTMemSize)
//                                     //
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  ===== Targets (TCx, xGA) for APOBEC mutations  =====================
    long TCxMemSize;    // = sizeof (TCXtag)
    char *TCXtag;       //contains tags from '*Xtag' with TCx_TAG
                        //                      for current segment
    long TCXsize;      // = real no Targets in *TCXtag  (<=TCxMemSize)
//
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  =====  Will be swiched before RNDcycles between
    //          XRObody (run without RT);  and RTsegment (RT run) =======
    char *pTagBody;     // may point to: a)Xtag or  b)RTtag
    long TagSize;       // may be eql :  a)Xsize or b)RTsize
//
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  GENES ================
    vector <HUGEN> genes;   // may be used for 'genes' or 'vecMuZon'
                            // but NOT TOGETHER !!!!!!!!!!!!!!!!!!!!

    vector < MUZONE > vecMuZon;
//
    long nNucInZon;  // n.Nucleas in Mut.Zones = size Mutation_Zones at curRTsegm
    long nNucOutCurZon; // n.Nucls outside of only current  'ZoneID'
    long nNucOutAllZon; // n.Nucls outside of ALL 'ZoneID'
//
    long nTCxInZon;      // n.Targets(TCx) in Mut.Zones of curRTsegm
    long nTCxOutCurZon;  //    "       "   outside of only current  'ZoneID'
    long nTCxOutAllZon;  //     "       "  outside of ALL 'ZoneID'

///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};  //curXpos=NULL;
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};    //curXpos=NULL;
    
	long read_DNK( const char *pFolder, char *fName );
    long markTCX_TAG ( );
    long read_Xzones( const char *pFolder, char *fName, int TagOnly=0); //  TagOnly=1 :: SET_MUZTWO_TAG
                                            //  not create 'vecMuZon'
        int testValidDNK( int Pos, char chREF, char chALT );
    int setZoneXtag( ); //PrepareBV( );
    int APOtest( long Pos, char chREF, char chALT );
    int creatTCx( char *pSelBuff );
    int creatRTSEGbody( long segStart, long segStop, char *pSelBuff );
    void clearRTSEG_TAG ( int noSeg );
    void markRTtag ( int setRT,  long First, long Last );
    void clcUnZonesArea( );
};
/////////////////////////////////////////////////////////

#define MUZONE_TAG  0x80
#define TCx_TAG     0x40
#define MUT_TAG     0x20        // real mutation tag

#define MUZTWO_TAG  0x08        // for save all zone tags together
#define RTSEG_TAG   0x01

#define SET_MUZONE_TAG(p_TAG)    *(p_TAG) |= MUZONE_TAG
#define GET_MUZONE_TAG(p_TAG)  (((*(p_TAG) & MUZONE_TAG)==0) ? 0 : 1 )

#define SET_MUZTWO_TAG(p_TAG)    *(p_TAG) |= MUZTWO_TAG
#define GET_MUZTWO_TAG(p_TAG)  (((*(p_TAG) & MUZTWO_TAG)==0) ? 0 : 1 )

#define SET_TCx_TAG(p_TAG)    *(p_TAG) |= TCx_TAG
#define GET_TCx_TAG(p_TAG)  (((*(p_TAG) & TCx_TAG)==0) ? 0 : 1 )

#define SET_MUT_TAG(p_TAG)    *(p_TAG) |= MUT_TAG
#define GET_MUT_TAG(p_TAG)  (((*(p_TAG) & MUT_TAG)==0) ? 0 : 1 )

#define SET_RTSEG_TAG(p_TAG)    *(p_TAG) |= RTSEG_TAG
#define GET_RTSEG_TAG(p_TAG)  (((*(p_TAG) & RTSEG_TAG)==0) ? 0 : 1 )

void  switchToXRO( );    // pointXRO( ); //
void  switchToRT( );     // pointRT( );  //

///////////////////////////////////////////////

int findXroByID( int xID );
int getNucID(const char Nuc);

#endif
