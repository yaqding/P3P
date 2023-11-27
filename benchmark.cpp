

#include "benchmark.h"

#include "lambdatwist/lambdatwist.p3p.h"
#include "new/new.p3p.h"
#include "new/new.p3p_direct.h"
#include "new/new.p3p_method1.h"
#include "new/new.p3p_method2.h"
#include "kneip/kneip.h"
#include "kneip/kneip2.h"
#include "utils/mlibtime.h"
#include <ke/ke.h>
// #include "/home/ding/lund/lambdatwist/jp/jp.h"
// #include "/home/ding/lund/lambdatwist/jp/jp_utils.h"
#include <p3p_generator.h>
#include <ke/ke_utils.h>
#include <data.h>
#include <kneip/kneip_utils.h>
#include <kneip/kneip_utils2.h>

#include "utils/string_helpers.h"

#include <fstream>

/*********
 *
 *
 *
 * Standardize input output format,
 * The speed comparison is primarily between two methods that output rotation matrixes. Therefore use that as basis.
 * The ke uses a silly format to maximize their speed...
 * convert to something resonable on output since it would always be followed by that anyways.
 *
 * Internally the methods should test the feasibility of the returned solutions, given the three points.
 * Ie dont return geometrically invalid solutions, nans, non rotations, transforms that move points behind cameras, or large reprojection errors
 * If they fail to do this it will always be done directly after.
 *
 * It would be very natural to use inheritance to create a nice modular structure for the benchmark. Doing so would hide the performance however,
 *  as vtable lookups take quite some time compared to the 300ns or so that the routines take to run... yes I tried it...
 *
 *
 *   Terms:
 * a valid solution is one considered valid by the solver, should always atleast test for geometric feasibility and nans
 * a correct solution is one which is both valid and satisfies the rotation matrix and reprojection criteria.
 *
 *
 * **/
#include <sstream>
#include <solver.h>

using std::cout;
using std::endl;
using namespace cvl;

namespace mlib
{
    namespace klas
    {

        template <class T, uint Rows, uint Cols>
        std::string displayMatrix(cvl::Matrix<T, Rows, Cols> M, bool displayrowcol)
        {

            std::vector<std::string> headers, rownames;
            if (displayrowcol)
            {
                for (uint i = 0; i < Cols; ++i)
                    headers.push_back(mlib::toStr(i));
                for (uint i = 0; i < Rows; ++i)
                    rownames.push_back(mlib::toStr(i) + ": ");
            }
            std::vector<std::vector<T>> rows;
            for (uint r = 0; r < Rows; ++r)
            {
                std::vector<T> row;
                for (uint c = 0; c < Cols; ++c)
                    row.push_back(M(r, c));
                rows.push_back(row);
            }

            return mlib::displayTable(headers, rows, rownames);
        }

        template <class T>
        class KeSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return kes::ke_p3p_fair(data, Rs, Ts);
            }
            std::string get_name() { return "Ke"; }
        };
        template <class T>
        class KneipSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return kneip::kneip_p3p_fair(data, Rs, Ts);
            }
            std::string get_name() { return "Kneip"; }
        };

        template <class T>
        class KneipSolver2
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {
                return kneip2::kneip_p3p_fair2(data, Rs, Ts);
            }
            std::string get_name() { return "nakano"; }
        };

        template <class T>
        class LambdaSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_lambdatwist(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "Lambda"; }
        };

        template <class T>
        class NewSolver
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "New"; }
        };

        template <class T>
        class NewSolverd
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "direct"; }
        };

        template <class T>
        class NewSolver1
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_method1(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "method1"; }
        };

        template <class T>
        class NewSolver2
        {
        public:
            int solve(Data<T> &data, cvl::Vector<Matrix<T, 3, 3>, 4> &Rs,
                      cvl::Vector<Vector<T, 3>, 4> &Ts)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                return p3p_new_method2(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            std::string get_name() { return "method2"; }
        };

        // verification of answers can be made slowly, so...

        template <class T, class Solver>
        P3PResult compute_accuracy(Solver S,
                                   std::vector<Data<T>> datas,
                                   T error_limit = 1e-6)
        {
            P3PResult res(S.get_name(), datas.size());

            std::ofstream outFile;
            //打开文件
            // outFile.open("/home/ding/Documents/Test.txt",std::ios::app);
            // outFile << 11111 <<"\n";

            for (Data<T> &data : datas)
            {
                cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                cvl::Vector<Vector<T, 3>, 4> Ts;

                int valid = S.solve(data, Rs, Ts);
                res.valid += valid; // valid according to solver.

                int duplicates = 0;
                int sols = data.good_solutions(Rs, Ts, valid, duplicates); // correct, with duplicates included

                if (sols == 0)
                    res.no_solution++;

                // if(sols==0)
                // {

                //     Vector3<Vector3<T>> yss=data.xr;
                //     Vector3<Vector3<T>> xss=data.x0;

                //     Pose<T> Pt = data.P;
                //     Vector3<T> tgt=Pt.t;
                //     Vector3<T> yss0 = yss[0]; Vector3<T> yss1 = yss[1]; Vector3<T> yss2 = yss[2];
                //     Vector3<T> xss0 = xss[0]; Vector3<T> xss1 = xss[1]; Vector3<T> xss2 = xss[2];
                //     cout.precision(17);
                //     // cout<<"error:"<<error<<endl;
                //     cout<<"data:"<<yss0(0)<<" "<<yss0(1)<<" "<<yss0(2)<<" "<<xss0(0)<<" "<<xss0(1)<<" "<<xss0(2)<<endl;

                //     cout<<"data:"<<yss1(0)<<" "<<yss1(1)<<" "<<yss1(2)<<" "<<xss1(0)<<" "<<xss1(1)<<" "<<xss1(2)<<endl;

                //     cout<<"data:"<<yss2(0)<<" "<<yss2(1)<<" "<<yss2(2)<<" "<<xss2(0)<<" "<<xss2(1)<<" "<<xss2(2)<<endl;
                //     cout<<"sols:"<<sols<<endl;
                //     cout<<"gt:"<<tgt(0)<<" " << tgt(1)<<" " << tgt(2)<<endl;

                //     for(int ii=0; ii<valid; ii++)
                //     {
                //         Vector3<T> Ts1=Ts[ii];
                //         cout<<"T:"<<Ts1(0)<<" "<<Ts1(1)<<" "<<Ts1(2)<<endl;
                //     }

                // }

                if (valid > sols)
                    res.incorrect_valid += (valid - sols);

                res.solutions += (sols - duplicates); // correct unique
                res.duplicates += duplicates;

                T error = data.min_error(Rs, Ts, valid);

                // outFile << error <<"\n";

                res.errors.push_back(error);
                if (error < error_limit)
                    res.ground_truth_in_set++;
                // else
                // {
                // Vector3<Vector3<T>> yss=data.xr;
                // Vector3<Vector3<T>> xss=data.x0;

                // Pose<T> Pt = data.P;
                // Vector3<T> tgt=Pt.t;
                //     Vector3<T> yss0 = yss[0]; Vector3<T> yss1 = yss[1]; Vector3<T> yss2 = yss[2];
                //     Vector3<T> xss0 = xss[0]; Vector3<T> xss1 = xss[1]; Vector3<T> xss2 = xss[2];
                //     cout.precision(17);
                //     cout<<"error:"<<error<<endl;
                //     cout<<"data:"<<yss0(0)<<" "<<yss0(1)<<" "<<yss0(2)<<" "<<xss0(0)<<" "<<xss0(1)<<" "<<xss0(2)<<endl;

                //     cout<<"data:"<<yss1(0)<<" "<<yss1(1)<<" "<<yss1(2)<<" "<<xss1(0)<<" "<<xss1(1)<<" "<<xss1(2)<<endl;

                //     cout<<"data:"<<yss2(0)<<" "<<yss2(1)<<" "<<yss2(2)<<" "<<xss2(0)<<" "<<xss2(1)<<" "<<xss2(2)<<endl;
                //     cout<<"sols:"<<sols<<endl;
                // cout<<"gt:"<<tgt(0)<<" " << tgt(1)<<" " << tgt(2)<<endl;

                //     for(int ii=0; ii<valid; ii++)
                //     {
                //         Vector3<T> Ts1=Ts[ii];
                //         cout<<"T:"<<Ts1(0)<<" "<<Ts1(1)<<" "<<Ts1(2)<<endl;
                //     }

                // }
            }
            outFile.close();
            return res;
        }

        template <class T>
        void test(const std::vector<Data<T>> &datas)
        {

            cout << "Beginning Test: " << endl;
            T error_limit = 1e-6;

            // P3PResult newp3p = compute_accuracy<T,NewSolver<T>>(NewSolver<T>(),datas,error_limit);
            P3PResult newp3p_direct = compute_accuracy<T, NewSolverd<T>>(NewSolverd<T>(), datas, error_limit);
            P3PResult newp3p_method1 = compute_accuracy<T, NewSolver1<T>>(NewSolver1<T>(), datas, error_limit);
            P3PResult newp3p_method2 = compute_accuracy<T, NewSolver2<T>>(NewSolver2<T>(), datas, error_limit);
            P3PResult lambda = compute_accuracy<T,LambdaSolver<T>>(LambdaSolver<T>(),datas,error_limit);
            P3PResult ke = compute_accuracy<T,KeSolver<T>>(KeSolver<T>(),datas,error_limit);
            P3PResult kneip= compute_accuracy<T,KneipSolver<T>>(KneipSolver<T>(),datas,error_limit);
            P3PResult kneip2= compute_accuracy<T,KneipSolver2<T>>(KneipSolver2<T>(),datas,error_limit);

            // make the result table
            std::vector<P3PResult> res = {newp3p_direct, newp3p_method1, newp3p_method2,lambda,ke,kneip,kneip2};
            // std::vector<P3PResult> res={newp3p,lambda,ke,kneip,kneip2};
            std::vector<std::string> headers;
            std::vector<std::string> columns;
            std::vector<std::vector<std::string>> rows;

            for (auto &r : res)
                headers.push_back(r.name);

            {
                // number of times ground truth was found
                std::vector<int> gt;
                for (auto &r : res)
                    gt.push_back(r.ground_truth_in_set);
                rows.push_back(toStrVec(gt));
                columns.push_back("ground truth found");
            }
            {
                // number of times a valid solution was found
                std::vector<int> gt;
                for (auto &r : res)
                    gt.push_back(datas.size() - r.no_solution);
                rows.push_back(toStrVec(gt));
                columns.push_back("any solution found");
            }
            {
                // ratio for ground truth found
                std::vector<double> gtratio;
                for (auto &r : res)
                    gtratio.push_back((double)r.ground_truth_in_set / datas.size());
                rows.push_back(toStrVec(gtratio));
                columns.push_back("ground truth found ratio");
            }

            {
                // number of no solution at all found
                std::vector<int> nosol;
                for (auto &r : res)
                    nosol.push_back(r.no_solution);
                rows.push_back(toStrVec(nosol));
                columns.push_back("no solution found");
            }
            {
                // number of valid solutions according to solver
                std::vector<int> solutions;
                for (auto &r : res)
                    solutions.push_back(r.valid);
                rows.push_back(toStrVec(solutions));
                columns.push_back("valid according to solver");
            }
            {
                // duplicates
                std::vector<int> duplicates;
                for (auto &r : res)
                    duplicates.push_back(r.duplicates);
                rows.push_back(toStrVec(duplicates));
                columns.push_back("duplicates");
                //
            }
            {
                // unique correct
                std::vector<int> usols;
                for (auto &r : res)
                    usols.push_back(r.solutions);
                rows.push_back(toStrVec(usols));
                columns.push_back("unique correct solutions");
                //
            }
            {
                // incorrect valid solutions
                std::vector<int> solutions;
                for (auto &r : res)
                    solutions.push_back(r.incorrect_valid);
                rows.push_back(toStrVec(solutions));
                columns.push_back("incorrect solutuons output by the solver");
            }

            cout << displayTable(headers, rows, columns, "Table: ") << endl;
            {
                std::vector<std::vector<float128>> errors;
                for (P3PResult r : res)
                    errors.push_back(r.errors);
                std::string str = display(errors, headers);
            }

            //         cout << "\033[2J\033[1;1H"; cout.flush();
        }

        template <class T>
        void versus(const std::vector<Data<T>> &datas,
                    int repeat_versus = 2,
                    bool inner_timers = false)
        {

            cout << "Beginning Versus" << endl;

            Timer kneip("Kneip: ", datas.size());
            Timer kneip2("nakano: ", datas.size());
            Timer ns("newsolver: ", datas.size());
            Timer lt("Lambda: ", datas.size());
            Timer kes("ke ", datas.size());
            Timer nsd("direct: ", datas.size());
            Timer ns1("m1: ", datas.size());
            Timer ns2("m2: ", datas.size());

            Timer tot_new("New Total");
            Timer tot_lambda("Lambda Total");
            Timer tot_kes("ke Total");
            Timer tot_kneip("Kneip Total");
            Timer tot_kneip2("nakano Total");
            Timer tot_newd("direct Total");
            Timer tot_new1("m1 Total");
            Timer tot_new2("m2 Total");

            // dont opt away
            int sols = 0;
            // ke

            bool dokneip = true;
            for (int i = 0; i < repeat_versus; ++i)
            {

                // kneip
                 if(dokneip){

                    tot_kneip.tic();
                    for(uint i=0;i<datas.size();++i){
                        //if(inner_timers){if(i==100){                kneip.clear();            }kneip.tic();}

                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+= kneip::kneip_p3p_fair(datas[i],Rs, Ts);

                        //    if(inner_timers)kneip.toc();
                    }

                    tot_kneip.toc();
                }
                // kneip2
                if(dokneip){

                    tot_kneip2.tic();
                    for(uint i=0;i<datas.size();++i){
                        //if(inner_timers){if(i==100){                kneip.clear();            }kneip.tic();}

                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+= kneip2::kneip_p3p_fair2(datas[i],Rs, Ts);

                        //    if(inner_timers)kneip.toc();
                    }

                    tot_kneip2.toc();
                }
                if(dokneip){
                    tot_kes.tic();

                    for(uint i=0;i<datas.size();++i){
                        //    if(inner_timers)if(i==100){            kes.clear();        }        kes.tic();
                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        sols+=kes::ke_p3p_fair(datas[i],Rs,Ts);
                        //      if(inner_timers)kes.toc();

                    }
                    tot_kes.toc();
                }
                // {
                //     tot_new.tic();
                //     for (uint i = 0; i < datas.size(); ++i)
                //     {
                //         //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                //         cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                //         cvl::Vector<Vector<T, 3>, 4> Ts;
                //         Data<T> data = datas[i];
                //         Vector3<Vector<T, 3>> yss = data.xr;
                //         Vector3<Vector<T, 3>> xss = data.x0;
                //         sols += p3p_new(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                //         //  if(inner_timers)lt.toc();
                //     }
                //     tot_new.toc();
                // }
                {
                    tot_newd.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_direct(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_newd.toc();
                }
                {
                    tot_new1.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_method1(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_new1.toc();
                }
                {
                    tot_new2.tic();
                    for (uint i = 0; i < datas.size(); ++i)
                    {
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
                        cvl::Vector<Vector<T, 3>, 4> Ts;
                        Data<T> data = datas[i];
                        Vector3<Vector<T, 3>> yss = data.xr;
                        Vector3<Vector<T, 3>> xss = data.x0;
                        sols += p3p_new_method2(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_new2.toc();
                }
                {
                    tot_lambda.tic();
                    for(uint i=0;i<datas.size();++i){
                        //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                        cvl::Vector<Matrix<T,3,3>,4> Rs;
                        cvl::Vector<Vector<T,3>,4> Ts;
                        Data<T> data=datas[i];
                        Vector3<Vector<T,3>> yss=data.xr;
                        Vector3<Vector<T,3>> xss=data.x0;
                        sols+=p3p_lambdatwist(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);
                        //  if(inner_timers)lt.toc();
                    }
                    tot_lambda.toc();
                }

                std::vector<mlib::Timer> tots = {tot_new, tot_lambda, tot_kes, tot_kneip, tot_kneip2,tot_newd,tot_new1,tot_new2};

                printtimers();
                std::vector<mlib::Timer> ts = {ns, lt, kes, kneip, kneip2,nsd,ns1,ns2};
                cout << tots << endl
                     << "wrote the ts" << endl;
            }
            if (false)
            {

                std::ofstream out("timing.csv");

                {
                    auto ts = kneip.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                {
                    auto ts = lt.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                {
                    auto ts = kes.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                {
                    auto ts = ns.getTimes();
                    for (Time time : ts)
                    {

                        out << (int)(time.ns) << ", ";
                    }
                    out << "\n";
                }
                out << std::endl;
                out.close();
                // cout<<"means: "<<mean(kn)<<" "<<mean(ls)<<endl;
            }

            // cout<<mlib::display(samples,false)<<endl;
            cout << "Versus Done !" << endl;
        }

        template <class T>
        void test_type()
        {
            // generate data
            cout << "generating";
            cout.flush();
            mlib::Timer timer("generator");

            Generator gennie;

            int N = std::pow(10, 7);

            std::vector<Data<T>> datas;
            datas.reserve(N);
            timer.tic();
            for (int i = 0; i < N; ++i)
            {

                datas.push_back(gennie.next<T>());

                // if((i%1000)==999)        cout<<"i: "<<i<<endl;
            }
            timer.toc();

            cout << timer << " done " << datas.size() << endl;

            test(datas);
            // printtimers();

            versus(datas);

            // testReal();
        }

        void test_special()
        {
            // generate data
            cout << "generating";
            cout.flush();
            mlib::Timer timer("generator");

            Generator gennie;

            int N = std::pow(10, 4);

            std::vector<Data<double>> datas;
            datas.reserve(N);
            timer.tic();
            for (int i = 0; i < N; ++i)
            {

                datas.push_back(gennie.special_case0<double>());

                // if((i%1000)==999)        cout<<"i: "<<i<<endl;
            }
            timer.toc();

            cout << timer << " done " << datas.size() << endl;

            test(datas);
            // printtimers();

            versus(datas);

            // testReal();
        }

        template <class T>
        void profile_lambda()
        {

            Generator gennie;

            int sols = 0;
            cvl::Vector<Matrix<T, 3, 3>, 4> Rs;
            cvl::Vector<Vector<T, 3>, 4> Ts;
            Data<T> data = gennie.next<T>();

            for (int i = 0; i < 100000; ++i)
            {

                Vector3<Vector3<T>> yss = data.xr;
                Vector3<Vector3<T>> xss = data.x0;
                sols += p3p_lambdatwist(yss[0], yss[1], yss[2], xss[0], xss[1], xss[2], Rs, Ts);
            }
            cout << sols << endl;
        }

        void testAll()
        {

            test_type<double>();
            // test_special();
            // profile_lambda<double>();
        }
    } // end namespace klas
} // end namespace mlib
