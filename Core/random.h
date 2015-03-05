//
//  random.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 1/29/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_random_h
#define Metaphor_random_h

namespace im
{
    // Simple random class that wraps common random operations
    class Rand
    {
    public:
        Rand() { time_seed(); }
        Rand(int seed_val) { seed(seed_val); }
        
        void time_seed() { m_rnd.seed((unsigned int)(std::chrono::system_clock::now().time_since_epoch().count())); }
        
        void seed(int seed_val)
        {
            if(seed_val<0)
                time_seed();
            else
                m_rnd.seed(seed_val==0 ? 1 : seed_val);
        }
        
        // includes both lower and upper in case of int
        int uniform_int(int lower, int upper)
        {
            std::uniform_int_distribution<int> uniform_dist(lower, upper);
            return uniform_dist(m_rnd);
        }
        
        // range is from lower to upper, not inluding upper
        float uniform_real(float lower, float upper)
        {
            std::uniform_real_distribution<float> uniform_dist(lower, upper);
            return uniform_dist(m_rnd);
        }
        
        // range is from lower to upper, not inluding upper
        double uniform_real(double lower, double upper)
        {
            std::uniform_real_distribution<double> uniform_dist(lower, upper);
            return uniform_dist(m_rnd);
        }
        
        // picks a random point within a rectangle
        MtxLoc point_in_rect(MtxRect const &rct)
        {
            return MtxLoc(uniform_int(rct.origin.row, rct.origin.row + rct.size.rows-1),
                            uniform_int(rct.origin.col, rct.origin.col + rct.size.cols-1));
        }

        // generate a normally distributed random variable with mean = 0 and std dev = 1
        double gauss() { return m_normal(m_rnd); }
        
    private:
        std::mt19937 m_rnd; // mersenne twister
        
        std::normal_distribution<double> m_normal;
    };
}

#endif
