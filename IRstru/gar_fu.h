//
//  gar_fu.h
//  garvard
//
//  Created by Ирина Пономарева on 05/04/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#ifndef gar_fu_h
#define gar_fu_h

/*
#define _C 0
#define _G 1
#define _A 2
#define _T 3

#define _C_N 0
#define _C_G 1
#define _C_A 2
#define _C_T 3

#define _C_N 0
*/
#define SPACER 3
#define START_PIN_SIZE 6

struct G_PIN {
    int length;
    int cnts[4][4];
    G_PIN() { length=0;
        for ( int r=0; r<4; r++)
            for ( int c=0; c<4; c++) cnts[r][c]=0;
    }
    G_PIN(int L) { length=L;
            for ( int r=0; r<4; r++)
                for ( int c=0; c<4; c++) cnts[r][c]=0;
    };
};


//////////////////////////////////////////////////////////////////

bool leseqv_PIN ( const G_PIN &x1, const G_PIN &x2 );

int defPrioNuc(char *pXb);
int defPinLen(char *pXb, char *First);

void printIRCap( FILE *fLoopRez, FILE *fRez );
void printIRsampl( FILE *fLoopRez, FILE *fRez, int nSamp );

#endif /* gar_fu_h */
