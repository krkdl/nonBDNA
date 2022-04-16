//
//  apobec.cpp
//  mutas
//
//  Created by Ирина Пономарева on 12/01/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mutas.h"
#include "xrosoma.h"

#ifdef DEBUG_TRACE_TSLD
extern FILE *tsld;
#endif

extern vector < XROSOMA > vecDNK;
extern vector < SAMPLE >  vecSAMPL;

//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int XROSOMA:: GQ_APOtest( long Pos, char chREF, char chALT )
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

