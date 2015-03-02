//
//  fft.h
//  ImpressiveAI
//
//  Created by SIMON WINDER on 2/26/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#ifndef ImpressiveAI_fft_h
#define ImpressiveAI_fft_h

namespace im
{
    enum FFTDirection
    {
        FFTDirectionForward = -1, FFTDirectionInverse = 1
    };
    
    template <typename TT>
    class FFTComplex
    {
    public:
        FFTComplex() : m_n(0) {}
        
        void init(int n);
    
        // Note that these functions do not scale the vector to normalize the transform
        void forward(VecView<std::complex<TT>> vv) { transform(vv, FFTDirectionForward); }
        void inverse(VecView<std::complex<TT>> vv) { transform(vv, FFTDirectionInverse); }
        void transform(VecView<std::complex<TT>> vv, FFTDirection direction);
        
    private:
        void factor_reduce(int &n, int factor);
        void transform_radix_2(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase);
        void transform_radix_3(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase);
        void transform_radix_4(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase);
        void transform_radix_5(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase);
        void transform_radix_N(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> vsrc, int product, int twiddlebase, int factor);
        void transform_power_of_2(VecView<std::complex<TT>> vv);
        
        int m_n;
        FFTDirection m_direction;
        std::vector<std::complex<TT>> m_twiddle;
        std::vector<std::complex<TT>> m_buffer;
        std::vector<int> m_factors;
        std::vector<int> m_offsets;
    };
    
    // Stand-alone FFT calls - these calls scale the vector for the inverse transform to make the composite trasform an identity
    
    // 1D:
    template <typename TT> void core_fft_forward(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc);
    template <typename TT> void core_fft_inverse(VecView<std::complex<TT>> vdst, VecView<std::complex<TT>> const &vsrc);
    template <typename TT> void core_fft_in_place_forward(VecView<std::complex<TT>> v);
    template <typename TT> void core_fft_in_place_inverse(VecView<std::complex<TT>> v);
    
    // 2D:
    template <typename TT> void core_fft_forward(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc);
    template <typename TT> void core_fft_inverse(MtxView<std::complex<TT>> mdst, MtxView<std::complex<TT>> const &msrc);
    template <typename TT> void core_fft_in_place_forward(MtxView<std::complex<TT>> m);
    template <typename TT> void core_fft_in_place_inverse(MtxView<std::complex<TT>> m);
    
    // Perform 1D FFTs of real-only data using efficient packing for faster transform.
    // Scaling is on inverse transform so that forward followed by inverse transform is the identity.
    // The size of the real data must be even and equal to twice the size of the complex FFT data vector.
    // The complex FFT vector includes three parts: the dc term packed in the real part of dst[0]. The N/2 frequency term packed in the imaginary part
    // of dst[0] and the positive frequency complex spectrum in the rest of the vector.
    // Source and destination are 1D vectors with n rows.
    template <typename TT> void core_fft_forward_real(VecView<std::complex<TT>> vdst, VecView<TT> const &vsrc);
    template <typename TT> void core_fft_inverse_real(VecView<TT> vdst, VecView<std::complex<TT>> const &vsrc);
}


#endif
