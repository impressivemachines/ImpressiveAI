//
//  mem.h
//  Metaphor Library
//
//  Created by SIMON WINDER on 2/4/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef Metaphor_mem_h
#define Metaphor_mem_h

namespace im
{
    // This class encapsulates a memory block which is aligned to the specified amount
    
    class MemoryBlock
    {
    public:
        MemoryBlock() : m_size(0), m_p(NULL) {}
        virtual ~MemoryBlock() { deallocate(); }
        
        // align_bytes must be a power of two
        uint8_t *resize(size_t size, size_t align_bytes = sizeof(void *))
        {
            // ensure that align_bytes is reasonable
            if(align_bytes < sizeof(void *))
                align_bytes = sizeof(void *);
            
            // ensure that size is a multiple of alignment
            size = align_bytes*((size + align_bytes - 1)/align_bytes);
            
            if(size==0)
                deallocate();
            else if(m_p==NULL || size!=m_size || !IM_IS_ALIGNED_K(m_p, align_bytes))
            {
                uint8_t *palloc;
                if(posix_memalign((void **)&palloc, align_bytes, size))
                    IM_THROW_ALLOC;
                
                if(m_p)
                    free(m_p);
                
                m_size = size;
                m_p = palloc;
            }
            
            return m_p;
        }
        
        void deallocate()
        {
            if(m_p)
                free(m_p);
            m_p = NULL;
            m_size = 0;
        }
        
        size_t size() { return m_size; }
        uint8_t const *ptr() const { return m_p; }
        uint8_t *ptr() { return m_p; }
        
    protected:
        size_t m_size;
        uint8_t *m_p;
    };
}

#endif
