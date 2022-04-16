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
//    int chrNum;
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
};
bool lesser_ZON ( const MUZONE &x1, const MUZONE &x2 );


//////////////////////////////////////////////////////////

class SAMPLE {
public:
		string SampName;
        vector < MUTANT > vMutNuc[NO_HUMAN_XRO]; // common set mutarions in XRO's

        int AcntSamp[NO_HUMAN_XRO];
        int BcntSamp[NO_HUMAN_XRO];
        SAMPLE(string S) { SampName=S;
            for (int n=0; n<NO_HUMAN_XRO; n++) { AcntSamp[n]=0; BcntSamp[n]=0; }
        };

        void calc_AB_APO( );
  //      void correctREAL_Xmut();
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
//  MUTATIONS ============================
    vector < MUZONE > vecMuZon;
//
//  AB_apobec ============================//
    int AcntTrg;
    int BcntTrg;
///////////
    XROSOMA() { chrNum=0; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};  //curXpos=NULL;
    XROSOMA(int cN) { chrNum=cN; Xsize=0L; Xbody=NULL; RTMemSize=0; TCxMemSize=0;
                TCXsize=0L; TCXtag=NULL; RTsize=0L;RTtag=NULL;};    //curXpos=NULL;
    
	long read_DNK( const char *pFolder, char *fName );
    long markTCX_TAG ( );
    long read_Xzones( const char *pFolder, char *fName);
        int testValidDNK( int Pos, char chREF, char chALT );
        int testSeqNuc(MUZONE &mz, char *buff);
    
    int AB_APOtest( long Pos, char chREF, char chALT );


};

/////////////////////////////////////////////////////////

#define A_TARG_TAG 0x80      //  [T|C]TCx  ||  xGA[G|A]
#define B_TARG_TAG 0x40      //  [G|A]TCx  ||  xGA[C|T]

#define SET_A_TARG(p_TAG)    *(p_TAG) |= A_TARG_TAG
#define GET_A_TARG(p_TAG)  (((*(p_TAG) & A_TARG_TAG)==0) ? 0 : 1 )

#define SET_B_TARG(p_TAG)    *(p_TAG) |= B_TARG_TAG
#define GET_B_TARG(p_TAG)  (((*(p_TAG) & B_TARG_TAG)==0) ? 0 : 1 )


int findXroByID( int xID );
int getNucIDn(const char Nuc);
int getCompNucID ( const char Nuc ); 


#endif
