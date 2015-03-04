//
//  geometry.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 3/1/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_geometry_h
#define ImpressiveAI_geometry_h

namespace im
{
    // Returns the coefficients of a line equation that passes through the thow points given
    // The coefficients vector has the form (a b c), for the line equation ax + by + c = 0
    template <typename TT> void core_line_through_points(VecView<TT> vcoefs, TT x1, TT y1, TT x2, TT y2);
}

#endif
