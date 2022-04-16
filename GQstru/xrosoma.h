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
    char Strand;
    MUZONE() { startPos=0; stopPos=0; leng=0; repeat=0; spacer=0; perm=0; Strand=0; };
};
bool lesser_ZON ( const MUZONE &x1, const MUZONE &x2 );

struct RTS_CNTS {   //   tCx  Strand  route     xGa  Strand  route
    int Stack_FW1;  //          +       1               -       0
    int Stack_BK0;  //          +       0               -       1
    int Stack_uDef; //          +       ?               -       ?
    int Line_FW1;   //          -       1               +       0
    int Line_BK0;   //          -       0               +       1
    int Line_uDef;  //          -       ?               +       ?
    RTS_CNTS() { Stack_FW1=0; Stack_BK0=0; Stack_uDef=0;
                 Line_FW1=0; Line_BK0=0; Line_uDef=0; };
};

//////////////////////////////////////////////////////////

class SAMPLE {
public:
		string SampName;
        vector < MUTANT > vMutNuc[NO_HUMAN_XRO]; // common set mutarions in XRO's
        RTS_CNTS cntAPO[NO_HUMAN_XRO];   
        SAMPLE(string S) { SampName=S;
        for (int n=0; n<NO_HUMAN_XRO; n++) { //APOinG_stack[n]=0; APOinC_strnd[n]=0; }
            cntAPO[n].Stack_FW1=0; cntAPO[n].Stack_BK0=0; cntAPO[n].Stack_uDef=0;
            cntAPO[n].Line_FW1=0;  cntAPO[n].Line_BK0=0;  cntAPO[n].Line_uDef=0; };
        };


        void markTCXtag(int nXro);
        int  calcGQ_APO( );
        void correctREAL_Xmut();
		void makeRndCycle( int nS, const char *canc_id);
        void printSampLine( FILE *fRezint);
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

    RTS_CNTS cntTARG;
//
///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};  //curXpos=NULL;
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};    //curXpos=NULL;
    
	long read_DNK( const char *pFolder, char *fName );
    long markTCX_TAG ( );
    long markTCX_TAG ( MUZONE &mz );
    long read_Xzones( const char *pFolder, char *fName, int TagOnly=0); //  TagOnly=1 :: SET_MUZTWO_TAG
                                                                        //  not create 'vecMuZon'
        int testValidDNK( int Pos, char chREF, char chALT );
        int testSeqNuc(MUZONE &mz, char *buff);
    void testZoneCross(  );
    int setZoneXtag( ); //PrepareBV( );
    int GQ_APOtest( long Pos, char chREF, char chALT );
    int creatRTSEGbody( long segStart, long segStop, char *pSelBuff );
    long calcTCx_OutZ ( );
    void makePINset ( );//vector< MUTANT> &vMutSet );
    int defPinMutation ( vector< MUTANT> &vMutSet );
    void iniPINset ( );
    
    int  defCrossZ( int curZ, int *sizCrossZ );
    int getCrossZ( int curZ );
    int  testGQ_ggg(long StartPos, long StopPos);
    void procCrossZones(long StartPos, long StopPos,
                   vector< MUTANT> &vMutSet, char *pDone );
    int procCornerC ( long StartPos, long indCorn, //vector < MUZONE >:: iterator iterZ,
                                vector< MUTANT> &vMutSet, char *pDoneIt );
    int markCornerC( );
};

/////////////////////////////////////////////////////////


#define MUZONE_TAG  0x80        // ZONE for forward Strand = '+'
#define TCx_TAG     0x40        // Target for APOBEC    (TGx + xGA)

#define MUZMNS_TAG  0x20        // ZONE for backward Strand = '-'
#define MUZTWO_TAG  0x08        // for save all zone tags together

#define RTS_FW      0x10        // replicator moves FORWARD
#define RTS_BK      0x01        // replicator moves BACK

#define SET_RTS_FW1(p_TAG)    *(p_TAG) |= RTS_FW
#define GET_RTS_FW1(p_TAG)  (((*(p_TAG) & RTS_FW)==0) ? 0 : 1 )

#define SET_RTS_BK0(p_TAG)    *(p_TAG) |= RTS_BK
#define GET_RTS_BK0(p_TAG)  (((*(p_TAG) & RTS_BK)==0) ? 0 : 1 )

#define SET_MUZONE_TAG(p_TAG)    *(p_TAG) |= MUZONE_TAG
#define GET_MUZONE_TAG(p_TAG)  (((*(p_TAG) & MUZONE_TAG)==0) ? 0 : 1 )

#define SET_MUZMNS_TAG(p_TAG)    *(p_TAG) |= MUZMNS_TAG
#define GET_MUZMNS_TAG(p_TAG)  (((*(p_TAG) & MUZMNS_TAG)==0) ? 0 : 1 )

#define SET_TCx_TAG(p_TAG)    *(p_TAG) |= TCx_TAG
#define GET_TCx_TAG(p_TAG)  (((*(p_TAG) & TCx_TAG)==0) ? 0 : 1 )

#define SET_MUZTWO_TAG(p_TAG)    *(p_TAG) |= MUZTWO_TAG
#define GET_MUZTWO_TAG(p_TAG)  (((*(p_TAG) & MUZTWO_TAG)==0) ? 0 : 1 )

int findXroByID( int xID );
int getNucIDn(const char Nuc);
int getCompNucID ( const char Nuc );
int  calcRoutStrand( char Nuc, char Strand, int Route, RTS_CNTS &cnts );

#endif
