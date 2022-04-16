#ifndef _XROSOMA_H__
#define _XROSOMA_H__


#include "vector"
#include "string"
using namespace std;

#define _C 0
#define _G 1
#define _A 2
#define _T 3

struct PROP_NUC {
    int prev;
    int ref;
    int alt;
    void ini() { prev=-1; ref=-1; alt=-1; };
};

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
			int noMut;              // Number mutation in curr_segment REAL case
			int noREAL_IN;          //  "     "  in mutZones  "   "    REAL case
			int noRAND_IN;          //  "     "  in mutZones  "   "    RANDOM case
			int cntLessReal;        //  counts where (sumRAND_IN <= sumREAL_IN)
			MUTATION() {  noMut=0; noREAL_IN=0; noRAND_IN=0; cntLessReal=0; };
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
		MUTATION mutSeg[NO_HUMAN_XRO];
            int sumMut;             // SUMM by Xro mutation in curr_segment REAL case
            int sumREAL_IN;         //   "      "  mut_in_Zone    "       " REAL case
            int sumRAND_IN;         //   "      "  mut_in_Zone    "       " RANDOM case
            int sumCntLess;		    //   "        counts where (sumRAND_IN <= sumREAL_IN)
                                    //  for all RND steps
        MUTATION mutAPO[NO_HUMAN_XRO];
            int sumMut_APO;         //
            int sumREAL_IN_APO;
            int sumRAND_IN_APO;
            int sumCntLess_APO;
    
        SAMPLE(string S) { SampName=S;
            sumMut=0; sumREAL_IN = 0; sumRAND_IN = 0; sumCntLess = 0;
            sumMut_APO=0; sumREAL_IN_APO = 0; sumRAND_IN_APO = 0; sumCntLess_APO = 0;
        };
        void clearSAMPLstat() {              // call first at each RT circl
            sumMut = 0;     sumMut_APO = 0;
            sumREAL_IN = 0; sumREAL_IN_APO = 0;
            sumRAND_IN = 0; sumRAND_IN_APO = 0;
            sumCntLess  =0; sumCntLess_APO =0;
            for (int n=0; n<NO_HUMAN_XRO; n++) {
                mutSeg[n].noMut=0;
                mutSeg[n].noREAL_IN=0;
                mutSeg[n].noRAND_IN=0;
                mutSeg[n].cntLessReal=0;
                mutAPO[n].noMut=0;
                mutAPO[n].noREAL_IN=0;
                mutAPO[n].noRAND_IN=0;
                mutAPO[n].cntLessReal=0;
            };
        };
		void clearRANDstat() {              // call first within each RANDOM circl
            sumRAND_IN = 0; sumRAND_IN_APO = 0;
            for (int n=0; n<NO_HUMAN_XRO; n++)  {
                mutSeg[n].noRAND_IN=0;
                mutAPO[n].noRAND_IN=0;
            };
        };
        void markTCXtag(int nXro);
        void calcREALstat( );
        void correctREAL_Xmut();
		void makeRndCycle( int nS, const char *canc_id);
        void printSamplLine( FILE *fRezint, int APO );
        void markMutations ( int ON );
};

void printTableCap( FILE *fRez );


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
    long StartSegIndx;      // Bounds of Current RTsegment in 'Xtag' (from 0 ....)
    long StopSegIndx;       // (from 0 ....)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  ===== Current RT segment. Will be updated for each segment. ================
    long RTMemSize;     // = sizeof (RTtag)
    char *RTtag;        // contains tags from '*Xtag' with RTSEG_TAG
    long RTsize;        // = real Segment_Size  (<=RTMemSize)
    long noNucInMuZon;//   XzoneSize;  // no Nucleas in Mut.Zones = size Mutation_Zones at curRTsegm
//                                     //
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//  ===== Targets (TCx, xGA) for APOBEC mutations  =====================
    long TCxMemSize;    // = sizeof (TCXtag)
    char *TCXtag;       //contains tags from '*Xtag' with TCx_TAG
                        //                      for current segment
    long TCXsize;      // = real no Targets in *TCXtag  (<=TCxMemSize)
    long noTCxInMuZon;    //TCXzoneSize;// no Targets in Mut.Zones current RTseg = sizeIN
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
    vector < MUZONE > vecMuZon;
//
//  IR_STUDY ============================
    int cntLoopMU[4][4];      // [prio] [alt]=[_ G A T ]
    int cntOtherMU[4][4][4][2];  // [prio] [ref] [alt]=[C G A T]  + [inLOOP]=[1 0]
//
///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};  //curXpos=NULL;
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};    //curXpos=NULL;
    
	long read_DNK( const char *pFolder, char *fName );
    long markTCX_TAG ( );
    long read_Xzones( const char *pFolder, char *fName);
        int testValidDNK( int Pos, char chREF, char chALT );
    int setZoneXtag( ); //PrepareBV( );
    int APOtest( long Pos, char chREF, char chALT );
    int creatTCx( char *pSelBuff );
    int creatRTSEGbody( long segStart, long segStop, char *pSelBuff );
    void clearRTSEG_TAG ( int noSeg );
    void markRTtag ( int setRT,  long First, long Last );
    void makePINset ( );//vector< MUTANT> &vMutSet );
    int defPinMutation ( vector< MUTANT> &vMutSet );
    void iniPINset ( );
    
    void  calcIRstat( vector< MUTANT> &vMutSet );
    int  defCrossZ( int curZ, int *sizCrossZ );
    long procCloop( vector < MUZONE >:: iterator iterZ, vector< MUTANT> &vMutSet );
    void procCrossZones(long StartPos, long StopPos,
                   vector< MUTANT> &vMutSet, char *pDone );
    int procCornerC ( long StartPos, long indCorn, //vector < MUZONE >:: iterator iterZ,
                                vector< MUTANT> &vMutSet, char *pDoneIt );
    int markCornerC( );
};
/////////////////////////////////////////////////////////


//#define NUCID(p_XBODY)         ( *(p_XBODY) & 0x7F )
#define MUZONE_TAG  0x80
#define TCx_TAG     0x40
#define MUT_TAG     0x20        // real mutation tag
#define CORN_TAG    0x10        // corner_C within Loop
#define LOOP_TAG    0x08        // Loop_part of mut_zone
#define RTSEG_TAG   0x01

#define SET_MUZONE_TAG(p_TAG)    *(p_TAG) |= MUZONE_TAG
#define GET_MUZONE_TAG(p_TAG)  (((*(p_TAG) & MUZONE_TAG)==0) ? 0 : 1 )

#define SET_TCx_TAG(p_TAG)    *(p_TAG) |= TCx_TAG
#define GET_TCx_TAG(p_TAG)  (((*(p_TAG) & TCx_TAG)==0) ? 0 : 1 )

#define SET_MUT_TAG(p_TAG)    *(p_TAG) |= MUT_TAG
#define GET_MUT_TAG(p_TAG)  (((*(p_TAG) & MUT_TAG)==0) ? 0 : 1 )

#define SET_RTSEG_TAG(p_TAG)    *(p_TAG) |= RTSEG_TAG
#define GET_RTSEG_TAG(p_TAG)  (((*(p_TAG) & RTSEG_TAG )==0) ? 0 : 1 )

#define SET_CORN_TAG(p_TAG)    *(p_TAG) |= CORN_TAG
#define GET_CORN_TAG(p_TAG)  (((*(p_TAG) & CORN_TAG)==0) ? 0 : 1 )

#define SET_LOOP_TAG(p_TAG)    *(p_TAG) |= LOOP_TAG
#define GET_LOOP_TAG(p_TAG)  (((*(p_TAG) & LOOP_TAG)==0) ? 0 : 1 )

void  switchToXRO( );    // pointXRO( ); //
void  switchToRT( );     // pointRT( );  //

///////////////////////////////////////////////

int findXroByID( int xID );
int getNucIDn(const char Nuc);

#endif
