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

struct MUT_SIGNA {                    // for ex.: AACC |  AGTT
    int lng;                          //            4      4
    int offs;                         //           -2     -1
    char sigset[NO_SIGNA_ID][8];      //
};
//struct APO_SIGN {       // for ex.:  tC  |  Ga
//    int lng;            //            2      2
//    int offs;           //           -1      0
//    char sign[8];       //           "TC"   "GA"
//};

class HUGEN {
public:
//    int chrNum;
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
//		int indXRO[NO_HUMAN_XRO];
        vector < MUTANT > vMutNuc[NO_HUMAN_XRO]; // common set mutarions in XRO's
		MUTATION muXr[NO_HUMAN_XRO];  //mutSeg
        MUTATION sumXr;
//            int sumMut;       // SUMM by Xro mutation in curr_segment REAL case
//            int sumREAL_IN;   //   "      "  mut_in_Zone    "       " REAL case
//            int sumRAND_IN;   //   "      "  mut_in_Zone    "       " RANDOM case
//            int sumCntLess;	//   "        counts where (sumRAND_IN <= sumREAL_IN)
                                    //  for all RND steps
        MUTATION muXrAPO[NO_HUMAN_XRO];  //mutAPO
        MUTATION sumXrAPO;
//            int sumMut_APO;         //
//            int sumREAL_IN_APO;
//            int sumRAND_IN_APO;
//            int sumCntLess_APO;
    
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
//    char *curXpos;
//
//  ===== real body of XRO will be read from files =============================
    long Xsize;
    char *Xbody;
    char *Xtag;     // 0x01=RTsegment; 0x80=inZone; 0x40=TCx; 0x20=APOBEC;
                    // RTsegment defines XRO work_area for creating TCXtag, noMut, REAL_IN, ets.
//    long StartSegIndx;      // Bounds of Current RTsegment in 'Xtag' (from 0 ....)
//    long StopSegIndx;       // (from 0 ....)
    
//                 dead zone; don't look at it.
    long StartNCBI;   //  beginZONE from 0
    long EndNCBI;     //  endZONE   from 0.
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
//  MUTATIONS ============================
//	vector < pair < unsigned long, unsigned long > >  vecMuZon;
                    //first=startPos; second=endPos;
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
        TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL; StartNCBI=-1; EndNCBI=-1;};
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;StartNCBI=-1; EndNCBI=-1;};
    
//	void iniXpos() { curXpos =Xbody; };
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


//#define NUCID(p_XBODY)         ( *(p_XBODY) & 0x7F )
#define MUZONE_TAG  0x10
#define TCx_TAG     0x40
#define MUT_TAG     0x20        // real mutation tag

#define MUZTWO_TAG  0x08        // for save all zone tags together
//#define RTSEG_TAG   0x01

#define SET_MUZONE_TAG(p_TAG)    *(p_TAG) |= MUZONE_TAG
#define GET_MUZONE_TAG(p_TAG)  (((*(p_TAG) & MUZONE_TAG)==0) ? 0 : 1 )

#define SET_MUZTWO_TAG(p_TAG)    *(p_TAG) |= MUZTWO_TAG
#define GET_MUZTWO_TAG(p_TAG)  (((*(p_TAG) & MUZTWO_TAG)==0) ? 0 : 1 )

#define SET_TCx_TAG(p_TAG)    *(p_TAG) |= TCx_TAG
#define GET_TCx_TAG(p_TAG)  (((*(p_TAG) & TCx_TAG)==0) ? 0 : 1 )

#define SET_MUT_TAG(p_TAG)    *(p_TAG) |= MUT_TAG
#define GET_MUT_TAG(p_TAG)  (((*(p_TAG) & MUT_TAG)==0) ? 0 : 1 )

//#define SET_RTSEG_TAG(p_TAG)    *(p_TAG) |= RTSEG_TAG
//#define GET_RTSEG_TAG(p_TAG)  (((*(p_TAG) & RTSEG_TAG)==0) ? 0 : 1 )

void  switchToXRO( );    // pointXRO( ); //
void  switchToRT( );     // pointRT( );  //

//#define IS_POS_INZONE_XRO(pos)  ( ( *(Xbody+pos) & 0x80)==0) )? 0 : 1 )
//#define IS_POS_INZONE_TCX(pos)  ( ( *(TCXbody+pos) & 0x80)==0) )? 0 : 1 )

//char nucID( char *pXb );
//void setZoneTag( char *pXb );
//int getZoneTag( char *pXb );
///////////////////////////////////////////////

//bool lesser ( const HUGEN &x1, const HUGEN &x2 );
//int LoadGenes( const char *pFolder, const char *fName );

int findXroByID( int xID );
int getNucID(const char Nuc);
int isIt_TCx( char *pXro );

#endif
