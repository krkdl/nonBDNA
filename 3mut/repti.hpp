//
//  repti.hpp
//  mutas
//
//  Created by Ирина Пономарева on 27/01/2021.
//  Copyright © 2021 Ирина Пономарева. All rights reserved.
//

#ifndef repti_hpp
#define repti_hpp

#include "mutas.h"
#include "xrosoma.h"

struct REPLYTIMESET {
    int indXro;
    unsigned long startPos;
    unsigned long stopPos;
    float repT;
    REPLYTIMESET(int a, unsigned long b, unsigned long c, float d)
    {indXro=a; startPos=b; stopPos=c; repT=d;}; 
    bool operator< (const REPLYTIMESET &a) const {
        return repT < a.repT;
    };
};


unsigned long readRTset( int indCancID );

#endif /* repti_hpp */
