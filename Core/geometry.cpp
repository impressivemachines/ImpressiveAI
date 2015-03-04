//
//  geometry.cpp
//  ImpressiveAI
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "imp_core.h"

template <typename TT>
void im::core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2)
{
    IM_CHECK_VALID(vcoefs);
    IM_CHECK_VECTOR_SIZE(vcoefs, 3);
    
    TT a = y1-y2;
    TT b = x2-x1;
    TT c = x1*y2 - x2*y1;
    
    TT norm = std::sqrt(a*a+b*b);
    if(norm==(TT)0)
    {
        vcoefs(0) = 0;
        vcoefs(1) = -1;
        vcoefs(2) = y1;
    }
    else
    {
        vcoefs(0) = a / norm;
        vcoefs(1) = b / norm;
        vcoefs(2) = c / norm;
    }
}

#define INST(TT) template void im::core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2)
INST(float); INST(double);
#undef INST
