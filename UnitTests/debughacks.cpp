//
//  debughacks.cpp
//  Metaphor
//
//  Created by SIMON WINDER on 3/5/15.
//  Copyright (c) 2015 Impressive Machines LLC. All rights reserved.
//

#include "meta_core.h"
#include "meta_vision.h"

void test1()
{
    float data[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    std::vector<float> vec;
    vec.resize(12);
    memcpy(&(vec[0]), data, vec.size() * sizeof(float));
    
    im::MtxView<float> mv1(3, 4, 4, 1, data);
    //mv1.print();
    
    im::MtxView<float> mv2(vec); // 1D
    //mv2.print();
    
    im::MtxView<float> mv3 = mv1.t();
    // mv3.print();
    
    // printf("%g\n",im::core_block_reduce_add(mv3));
    
    im::core_block_sort(mv3, im::SortDirectionDescending);
    // mv3.print();
    
    im::core_block_clear_lower_tri(mv3, false);
    //  mv3.print();
    
    //  float d2[4] = {10, 7, 2, 0 };
    //  im::MtxView<float> mv4(4,1,1,0,d2);
    
    //  printf("%d\n", im::core_block_nearest(mv3, mv4));
    
    //  float foo[100];
    //  im::MtxView<float> mv5(4,4,4,1,foo);
    //  mv5 = 0;
    //  im::core_block_matrix_multiply_ATA(mv5, mv3.t());
    // mv5.print();
    
    //  im::core_block_matrix_multiply(mv5, mv3, mv3.t());
    //  mv5.print();
    
    im::Rand rnd;
    rnd.seed(56);
    //   for(int i=0; i<10; i++)
    //       printf("%g\n", rnd.uniform_real(0.0f,1.0f));
    
    
    // printf("%lu\n", sizeof(im::TypeProperties<uint8_t>::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties<short>::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties<int>::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties<float>::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties<double>::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties< std::complex<float> >::Accurate));
    //  printf("%lu\n", sizeof(im::TypeProperties< std::complex<double> >::Accurate));
    
    im::Mtx<float> mf(4,4);
    mf.random_gaussian(rnd);
    //    mf.print();
    
    mf.save_binary("test.imp");
    
    im::Mtx<float> mg(4,4);
    mg.load_binary("test.imp");
    
    //  mg.print();
    
    mg << "[ 1 2 3 ; 5 6 7; 0;1 ]";
    //  mg.print();
    //  mg.print_size();
    
    //im::MtxViewIterator<float> it = mg.view().iterator_begin();
    //im::MtxViewIterator<float> itend = mg.view().iterator_end();
    
    
    im::MtxViewIterator<float> it = mg.view();
    
    // /  for(; it.is_active(); ++it)
    //   {
    //       printf("%d %d -> %g\n", it.row, it.col, *it);
    //   }
    
    it.row = 1;
    it.col = 2;
    
    //   do {
    //       printf("%d %d -> %g\n", it.row, it.col, it.at(-1, 0));
    //   } while(it.next());
    
    im::Mtx<float> mm(5,5);
    im::MtxViewIterator<float> it2(mm.view());
    do {
        *it2 = it2.row * 5 + it2.col + 1;
    } while(it2.next());
    
    // mm.print();
    
    im::Permute p;
    p.resize(5);
    p.swap(1, 4);
    p.swap(1, 2);
    p.swap(2, 3);
    p.swap(3, 1);
    p.swap(4, 0);
    
    im::Mtx<float> mpe = p.matrix_P<float>();
    //  mpe.print();
    
    p.matrix_PX_in_place(mm.view());
    //  mm.print();
    p.matrix_PTX_in_place(mm.view());
    //  mm.print();
    
    im::MatrixDecompLU<float> lu(mf.view());
    
    im::Mtx<float> inv = lu.inverse();
    // inv.print();
    // (inv*mf).print();
    
    im::Mtx<float> mfmf = mf.t() * mf;
    
    im::Mtx<float> mfx(4,3),mfy(4,3);
    //Ax = y
    mfx.random_gaussian(rnd);
    mfy = mfmf * mfx;
    
    
    im::MatrixDecompLLT<float> ldl(mfmf.view());
    
    // (ldl.inverse()*mfmf).print();
    
    //   ldl.solve(mfy.view()).print();
    //   mfx.print();
    
    //   mf.print();
    
    im::MatrixDecompQR<float> qr(mf.view());
    
    im::Mtx<float> uu = (qr.inverse() * mf);
    im::Mtx<float> uu1(uu.rows(), uu.cols());
    uu1.set_identity();
    //  uu1.print();
    //   uu.print();
    //printf("%d\n", im::core_equal(uu, uu1, 0.000001));
    
    mf.print();
    for(int r=0;r<4; r++)
        for(int c=0;c<4;c++)
        {
            uu(r,c) = mf.sample_bilinear(0.35f*r, 0.5f*c) - mf.sample_bicubic(0.35f*r, 0.5f*c);
        }
    uu.print();
}

void test2()
{
    im::Rand rnd;
    im::Md m1(10,10), m2(10,10);
    m1.random_gaussian(rnd);
    m2 = m1 * m1.t();
    
    printf("start\n");
    im::MatrixDecompTridiag<double> td(m2.view());
    printf("end\n");
    
    m2.print();
    
    printf("coefs\n");
    td.householder_packed().print();
    td.householder_coefs().print();
    
    
    im::Md mQ = td.matrix_Q();
    im::Md mT = td.matrix_T();
    
    printf("Q:\n");
    mQ.print();
    printf("T:\n");
    mT.print();
    
    printf("pair:\n");
    m2.print();
    printf("%g\n", ((mQ * mT * mQ.t()) - m2).max_abs());
    
    
}

void test3()
{
    im::Mf min(500,500);
    im::Rand rnd;
    min.random_gaussian(rnd);
    
    im::Mf mk(9,9);
    mk.random_gaussian(rnd);
    //mk = 0;
    //mk(2,2) = 1;
    // mk.print();
    // min.print();
    
    im::Mf mout(250,125);
    
    im::core_convolve(mout.view(), mk.view(), min.view(), im::PadModeExtend, im::MtxLoc(mk.rows()/2,mk.cols()/2), im::MtxLoc(0,0), 0.0f, 4, 2, 64);
    
    // mout.print();
    printf("done\n");
    
    for(int i=0; i<mout.rows(); i++)
        for(int j=0; j<mout.cols(); j++)
        {
            float sum = 0;
            for(int kr = 0; kr<mk.rows(); kr++)
                for(int kc=0; kc<mk.cols(); kc++)
                {
                    int krin = kr + i*4 - mk.rows()/2;
                    int kcin = kc + j*2 - mk.cols()/2;
                    krin = im::core_clamp(krin, min.rows());
                    kcin = im::core_clamp(kcin, min.cols());
                    //if(krin>=0 && krin<min.rows() && kcin>=0 && kcin<min.cols())
                    //    sum += mk(mk.rows() - 1 - kr,mk.cols() - 1 - kc) * min(krin,kcin);
                    sum += mk(mk.rows() - 1 - kr,mk.cols() - 1 - kc) * min(krin,kcin);
                }
            mout(i,j) -= sum;
        }
    printf("done2\n");
    //  mout.print();
    printf("%g\n", mout.max_abs());
    
}

void test4()
{
    im::Rand rnd;
    im::Vf vin(100000);
    vin.random_gaussian(rnd);
    im::Vf vout(50000);
    im::Vf vk(1023);
    vk.random_gaussian(rnd);
    //vk.print();
    //vin.print();
    
    im::core_convolve(vout.view(), vk.view(), vin.view(), im::PadModeExtend, vk.rows()/2, 0, 0.0f, 2);
    
    printf("done\n");
    //vout.print();
    for(int i=0; i<vout.rows(); i++)
    {
        float sum = 0;
        for(int kr=0; kr<vk.rows(); kr++)
        {
            int krin = kr + i*2 - vk.rows()/2;
            krin = im::core_clamp(krin, vin.rows());
            // if(krin>=0 && krin<vin.rows())
            sum += vk(vk.rows()-1-kr) * vin(krin);
        }
        vout(i) -= sum;
    }
    //vout.print();
    
    printf("%g\n", vout.max_abs());
    
}

void test5()
{
    int n = 5;
    
    im::Mtx<double> m1(n,n), m2(n,n);
    im::Rand rnd;
    //rnd.seed(23423);
    
    m1.random_gaussian(rnd);
    //  m1.row(0).copy_from(m1.row(1));
    
    m2 = m1*m1.t();
    
    im::MatrixDecompEigenSymmetric<double> eig(m2.view());
    
    eig.eigenvalues().print();
    //eig.eigenvectors().print();
    
    // m2.print();
    //printf("%g\n", (eig.eigenvectors() * eig.eigenvalues().diag_matrix() * eig.eigenvectors().t() - m2).max_abs());
    
    im::Vec<double> vr(n), vv(n);
    vv.random_gaussian(rnd);
    
    vr = eig.solve(vv.view());
    
    // ((m2 * vr) - vv).print();
    
    im::Mtx<double> mr(n,4), mv(n,4);
    mv.random_gaussian(rnd);
    mr = eig.solve(mv.view());
    
    printf("%g\n", ((m2 * mr) - mv).max_abs());
    
    im::Mtx<double> msq = eig.square_root();
    
    (msq * msq - m2).print();
    (msq * eig.inverse_square_root()).print();
}

void test6()
{
    int m = 4;
    int n = 9;
    
    im::Mtx<double> m1(m,n);
    im::Rand rnd;
    m1.random_gaussian(rnd);
    
    im::MatrixDecompQR<double> qr;
    qr.compute(m1.view(), false);
    
    qr.matrix_Q_thin().print();
    qr.matrix_R_thin().print();
    
    printf("error = %g\n",(qr.matrix_Q_thin() * qr.matrix_R_thin() - m1).max_abs());
    
    qr.matrix_Q().print();
    qr.matrix_R().print();
    
    printf("error = %g\n",(qr.matrix_Q() * qr.matrix_R() - m1).max_abs());
    
}

void test7()
{
    int m = 20;
    int n = 5;
    im::Mtx<double> m1(m,n), m2(m,n);
    im::Rand rnd;
    m1.random_gaussian(rnd);
    //m1.col(3).copy_from(m1.col(2));
    
    m2.copy_from(m1);
    
    im::MatrixDecompSVD<double> svd;
    
    
    svd.compute(m1.view(), true);
    
    printf("%g\n", (svd.matrixU() * svd.vectorS().diag_matrix() * svd.matrixV().t() - m2).max_abs());
    
    //im::Mtx<double> mq(n,n);
    //im::Vec<double> vs(n);
    
    
    //m1.print();
    
    // svd.in_place_jacobi(m1.view(), mq.view(), vs.view());
    
    
    //m1.print();
    // mq.print();
    //vs.print();
    
    //printf("%g\n", (m1 * vs.diag_matrix() * mq.t() - m2).max_abs());
    
    im::Mtx<double> mx(n,4), my(m,4), mz(n,4);
    mx.random_gaussian(rnd);
    
    my = m1 * mx;
    mz = svd.solve(my.view());
    printf("%g\n", (mz - mx).max_abs());
    
    svd.pseudo_inverse().print();
    (svd.pseudo_inverse() * m1).print();
    
    im::Mtx<double> mA(2,2), mV(2,2), mU(2,2);
    im::Vec<double> vS(2);
    mA.random_gaussian(rnd);
    //mA.col(1).copy_from(mA.col(0));
    
    core_decomp_svd_2x2(mU.view(), vS.view(), mV.view(), mA.view());
    
    mA.print();
    mU.print();
    vS.print();
    mV.print();
    
    (mU * vS.diag_matrix() * mV.t()).print();
}

void test8()
{
    im::Rand rnd;
    im::Mtx<float> m(2,30);
    
    for(int i=0; i<m.cols(); i++)
    {
        m(0,i) = i*2;
        m(1,i) = 7 + 3*i+ rnd.gauss()*(rnd.uniform_real(0.0f, 1.0f)>0.7 ? 10 : 0.01);
        // m(2,i) = i > 5 ? 5 : 1;
    }
    
    m.print();
    FILE *fp = fopen("a.txt", "w");
    for(int i=0; i<m.cols(); i++)
        fprintf(fp,"%g %g\n", m(0,i), m(1,i));
    fclose(fp);
    im::Vec<float> v(3);
    
    core_line_fit_2d(v.view(), m.view());
    v.print();
    
    core_line_fit_2d_orthog_robust(v.view(), m.view(), rnd);
    v(0)= -v(0)/v(1);
    v(2) = -v(2)/v(1);
    v(1) = -1;
    v.print();
    
    im::MtxView<float> mu = m.block(0,0,2,10).view();
    
    im::Vec<float> lfo(2), lfu(2);
    core_line_fit_orthog(lfo.view(), lfu.view(), mu);
    lfo.print();
    lfu.print();
    
    m.print();
    
    for(int i=0; i<m.cols(); i++)
    {
        printf("%g ", i*lfu(1)/lfu(0) - lfo(0)*lfu(1)/lfu(0) + lfo(1));
    }
    printf("\n");
    
    //return;
    
    v.print();
    m.print();
    
    im::Mtx<float> mv(2,2);
    im::Vec<float> vl(2);
    core_stats_pca(mv.view(), vl.view(), mu);
    
    vl.sqrt().print();
    mv.print();
    
    im::Vec<float> vmean = m.block(0,0,2,10).reduce_rows_mean();
    //vmean.print();
    for(int i=0; i<10; i++)
    {
        m(0,i) -= vmean(0);
        m(1,i) -= vmean(1);
    }
    im::MatrixDecompSVD<float> svd(mu.t());
    //svd.matrixU().print();
    svd.vectorS().print();
    svd.matrixV().print();
    
    float k = 1/sqrt(10);
    (svd.vectorS() * k).print();
    
    
    im::Mtx<float> mls(50,4);
    mls.random_gaussian(rnd);
    
    im::Vec<float> vls(50), xls(4);
    vls.random_gaussian(rnd);
    xls.random_gaussian(rnd);
    vls = mls * xls;
    for(int i=0; i<vls.rows(); i++)
        vls(i) += rnd.gauss() * 0.001;
    
    core_least_squares_normal(xls.view(), mls.view(), vls.view());
    xls.view().print();
    
    core_least_squares_svd(xls.view(), mls.view(), vls.view());
    xls.view().print();
    
    im::Mtx<float> ata = mls.t() * mls;
    im::Vec<float> atv = mls.t() * vls;
    // ata.print();
    // atv.print();
    
    im::MatrixDecompLDLT<float> ldlt(ata.view());
    
    im::Vec<float> sol = ldlt.solve(atv.view());
    sol.print();
    
    (ata * sol - atv).print();
    (mls * sol - vls).print();
    
    svd.compute(mls.view());
    svd.matrixV().print();
    svd.vectorS().print();
    printf("%g\n", (mls - svd.matrixU() * svd.vectorS().diag_matrix() * svd.matrixV().t()).max_abs());
    
    im::Vec<float> vtest(50), vresult(4);
    vresult.random_gaussian(rnd);
    vtest = mls * vresult;
    printf("%g\n", (svd.solve(vtest.view()) - vresult).max_abs());
    
    sol = svd.solve(vls.view());
    (mls * sol - vls).print();
    
}

void test9()
{
    im::Mtx<float> mX = "[ 8 6 2; 7 4 9; 4 8 2]";
    mX.inverse().print();
    
    im::Mtx<float> m(3,3), m4(4,4);
    im::Vec<float> v(3);
    v(0)=1;
    v(1) = 0;
    v(2) = 0;
    
    core_make_3x3_skew_symmetric(m.view(), 0.01f, 0.0f, 0.0f);
    //m.print();
    core_make_3x3_rotation_about_x(m.view(), 0.01f);
    m.print();
    core_make_3x3_skew_symmetric(m.view(), 0.0f, 0.1f, 0.0f);
    //m.print();
    core_make_3x3_rotation_about_y(m.view(), 0.01f);
    m.print();
    core_make_3x3_skew_symmetric(m.view(), 0.0f, 0.0f, 0.1f);
    //m.print();
    core_make_3x3_rotation_about_z(m.view(), 0.01f);
    m.print();
    core_make_3x3_rotation_euler(m.view(), 0.01f, 0.0f, 0.0f);
    //m.print();
    core_make_3x3_rotation_axis_angle(m.view(), v.view(), 0.01f);
    //m.print();
    core_make_4x4_rotation_about_x(m4.view(), 0.01f);
    //m4.print();
    
    core_make_3x3_rotation_euler(m.view(), 0.01f, 0.0f, 0.0f);
    m.print();
    
    core_make_3x3_rotation_euler(m.view(), 0.00f, 0.01f, 0.0f);
    m.print();
    
    core_make_3x3_rotation_euler(m.view(), 0.00f, 0.0f, 0.01f);
    m.print();
    
}

/*void test10()
{
    im::GenericImgView iv;
    float data[256];
    im::ImgView<float> ivf(8,16,2,16,data);
    
    ivf = 0.0f;
    ivf.block(0,0,3,3).fill(1.0f);
    ivf.matrix_view(0) = 42.0f;
    ivf.print_pixel(1, 2);
    ivf.print_pixel(3, 4);

    im::Img<uint8_t> myimg(128,64);
    myimg(0,0) = 128;
    myimg.print_size();
    myimg.print_pixel(0, 0);
    
  //  myimg.wrap(16, 16, 1, 16, data);
    printf("done\n");
}
*/
void test11()
{
    im::Quat<float> q1,q2,q3;
    
    q1.from_euler_angles(1, -1, -0.99);
    float r,p,y;
    q1.to_euler_angles(p,y,r);
    printf("%g %g %g\n", p, y, r);
    im::Mtx<float> mrot(3,3);
    q1.to_rotation_matrix(mrot.view());
    mrot.print();
    core_make_3x3_rotation_euler(mrot.view(), p, y, r);
    mrot.print();
    q1.from_rotation_matrix(mrot.view());
    q1.to_euler_angles(p,y,r);
    printf("%g %g %g\n", p, y, r);
    
    q1.from_euler_angles(0.1, 0, 0);
    q2.from_euler_angles(0.4, 0, 0);
    
    q3.from_iterpolation(q1, q2, 0.5f);
    q3.to_euler_angles(p, y, r);
    printf("%g %g %g\n", p, y, r);
    
    float v1[3], v2[3];
    v1[0] = 1;
    v1[1] = 0;
    v1[2] = 0;
    v2[0] = 1;
    v2[1] = 0;
    v2[2] = -1;
    q3.from_vector_rotation(im::VecView<float>(3, 1, v1), im::VecView<float>(3,1,v2));
    q3.to_rotation_matrix(mrot.view());
    mrot.print();
    q3.to_euler_angles(p, y, r);
    printf("%g %g %g\n", p*CONST_180_PI, y*CONST_180_PI, r*CONST_180_PI);
    
    q1.from_euler_angles(0.2, 0, 0);
    q2.from_euler_angles(0.5, 0.3, 0);
    q3 = q2 * q1;
    q3.to_euler_angles(p, y, r);
    printf("%g %g %g\n", p, y, r);
    float u[3];
    float ang;
    q3.to_axis_angle(im::VecView<float>(3,1,u), ang);
    printf("%g %g %g %g\n", u[0],u[1],u[2],ang);
    q3.from_axis_angle(im::VecView<float>(3,1,u), ang);
    q3.to_euler_angles(p,y,r);
    printf("%g %g %g\n", p, y, r);
}

class FuncTest : public im::FuncEval1D<double>
{
public:
    double eval_fx(double x)
    {
     //   printf("ev %.16f\n",x);
        return log(std::abs(x)+1)+sin(x)-2;
    }
};

class FuncTest2 : public im::FuncEval1D<double>
{
public:
    double eval_fx(double x)
    {
        printf("ev %.16f\n",x);
        return sin(x);
    }
    
    double eval_dfx(double x)
    {
        printf("evd %.16f\n",x);
        return cos(x);
    }
};

void test12()
{
   // double r = im::core_root_search(func, NULL, 1.0, 8.0);
   // printf("%.10g\n",r);
    
    FuncTest ft;
    FuncTest2 ft2;
    
    double xa, xb, xc;
    im::core_line_min_bracket(xa, xb, xc, &ft, 12.71, 12.74);
    printf("%g %g %g\n", xa, xb, xc);
    
    double xmin, fxmin;
    im::core_line_min(xmin, fxmin, &ft, 14.0, 18.0);
    printf("%g %g\n", xmin, fxmin);
    
    im::core_line_min_using_derivs(xmin, fxmin, &ft2, 0.0, 5.0);
    printf("%g %g\n", xmin, fxmin);
    
}

class LMFunc : public im::LevenbergMarquardt<double>
{
    void eval_residual(im::Vec<double> &vresidual, im::Vec<double> const &vx)
    {
        for(int i=0; i<vx.rows(); i++)
            vresidual(i) = 1-cos(vx(i));
    }
    
    void eval_jacobian(im::Mtx<double> &mjacobian, im::Vec<double> &vx)
    {
        for(int i=0; i<vx.rows(); i++)
            for(int j=0; j<vx.rows(); j++)
            {
                if(i==j)
                    mjacobian(i,j) = sin(vx(i));
                else
                    mjacobian(i,j) = 0;
            }
    }
};

void test13()
{
    LMFunc lm;
    
    im::Rand rnd;
    im::Vec<double> vx(5);
    vx.random_gaussian(rnd);
    vx *= 10.0;
    vx.print();
    lm.init(vx);
    
    while(!lm.step())
    {
        lm.state().print();
        printf("err=%g derr=%g dx=%g lam=%g\n", lm.error(), lm.delta_error(), lm.delta_x(), lm.lambda());
    }
    
    (lm.state()/CONST_PI).print();
    if(lm.early_exit())
        printf("failed\n");
}

class SGFunc : public im::StochasticMin<double>
{
public:
    void setup(int d)
    {
        im::Mtx<double> ma(d,d);
        im::Rand rnd;
        ma.random_gaussian(rnd);
        m = ma.t() * ma;
        
        im::MatrixDecompEigenSymmetric<double> eig(m.view());
        m_cond = eig.eigenvalues()(0) / eig.eigenvalues()(d-1);
    }
    
    double eval_fx(im::Vec<double> const &vx)
    {
        return vx.dot_product((m * vx));
    }
    
    void eval_dfx(im::Vec<double> &vdfx, im::Vec<double> const &vx)
    {
        im::Vec<double> grad = 2.0 * (m*vx);
        printf("grad = ");
        grad.print();
        vdfx.copy_from(grad);
    }
    
    /*double eval_fx(im::Vec<double> const &vx)
    {
        double v = 1;
        for(int i=0; i<vx.rows(); i++)
            v *= cos(vx(i));
        return v;
    }
    
    void eval_dfx(im::Vec<double> &vdfx, im::Vec<double> const &vx)
    {
        for(int i=0; i<vx.rows(); i++)
        {
            double v = 1;
            for(int j=0; j<vx.rows(); j++)
                if(i==j)
                    v *= -sin(vx(j));
                else
                    v *= cos(vx(j));
            vdfx(i) = v;
        }
    }*/
    
    im::Mtx<double> m;
    double m_cond;
};

void test14()
{
    SGFunc sg;
 
    int d = 5;
    sg.setup(d);
    im::Rand rnd;
    im::Vec<double> vx(d);
    vx.random_gaussian(rnd);
    vx *= 10.0;
    vx.print();
    sg.init(vx);
    
    while(!sg.step())
    {
        sg.state().print();
        printf("%d r=%g m=%g fx=%g dfx=%g dx=%g\n", sg.iteration_count(), sg.rate(), sg.momentum(), sg.fx(), sg.delta_fx(), sg.delta_x());
    }
    
    sg.state().print();
    printf("condition number was %g\n", sg.m_cond);
}

class PFunc : public im::PowellMin<double>
{
public:
    void setup(int d)
    {
        im::Mtx<double> ma(d,d);
        im::Rand rnd;
        ma.random_gaussian(rnd);
        m = ma.t() * ma;
        
        im::MatrixDecompEigenSymmetric<double> eig(m.view());
        m_cond = eig.eigenvalues()(0) / eig.eigenvalues()(d-1);
    }
    
    double eval_fx(im::Vec<double> const &vx)
    {
       // printf("eval: ");
       // vx.print();
    
        return vx.dot_product((m * vx));
    }
    
    im::Mtx<double> m;
    double m_cond;
};

void test15()
{
    PFunc powell;
    im::PowellMinParams params;
    params.line_min_eps = 1e-4;
    
    int d = 5;
    powell.set_parameters(params);
    powell.setup(d);
    
    im::Rand rnd;
    im::Vec<double> vx(d);
    vx.random_gaussian(rnd);
    vx *= 10.0;
    vx.print();
    powell.init(vx);
    
    while(!powell.step())
    {
        powell.state().print();
        printf("%d fx=%g dfx=%g dx=%g\n", powell.iteration_count(), powell.fx(), powell.delta_fx(), powell.delta_x());
    }
    
    powell.state().print();
    printf("condition number was %g\n", powell.m_cond);
}

int main()
{
    test15();
    return 0;
}

