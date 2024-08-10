using Microsoft.Win32.SafeHandles;
using System;
using System.Linq.Expressions;
using System.Security.Cryptography;
using System;

using Double5Func = System.Func<double, double, double, double, double, double>;
using Double3Func = System.Func<double, double, double, double>;
using Double2Func = System.Func<double>;
using System.Data.SqlTypes;
using System.Xml.Schema;

namespace lab_05
{
    internal class Model
    {
        public double a1, b1, c1, m1;
        public double beta, x0, z0, f0;
        public double hx, hz, tau;
        public double F0, T0, S, l;
        public double alpha2, alpha3, alpha4;
        public Model(double a1, double b1, double c1, double m1, double f0, double beta,
            double x0, double z0, double hx, double hz, double F0, double T0, double alpha2,
            double alpha3, double alpha4)
        {
            this.a1 = a1;
            this.b1 = b1;
            this.c1 = c1;
            this.m1 = m1;
            this.beta = beta;
            this.x0 = x0;
            this.z0 = z0;
            this.f0 = f0;
            this.hx = hx;
            this.hz = hz;
            this.F0 = F0;
            this.T0 = T0;
            this.alpha2 = alpha2;
            this.alpha3 = alpha3;
            this.alpha4 = alpha4;
            this.tau = 1;
            S = 30;
            l = 10;
        }

        public double Labmda(double T)
        {
            return 1.5;
            //return a1 * (b1 + c1 * Math.Pow(T, m1));
        }

        public double F(double x, double z)
        {
            return 0;
            // return f0 * Math.Exp(-beta * ((x - x0) * (x - x0) * (z - z0) * (z - z0)));
            //return f0 * Math.Exp(-beta * ((x - 3) * (x - 3) + (z - 3) * (z - 3)));
        }

        private double An()
        {
            double ln_12m = Labmda(1);
            return ln_12m / (hx * hx);
        }

        private double Bn()
        {
            return An() + Cn() + 2 / tau;
        }

        private double Cn()
        {
            double ln12m = Labmda(1);
            return ln12m / (hx * hx);
        }

        public double Dn(double ynm_1, double ynm, double ynm1, double x, double z)
        {
            double lnm12 = (Labmda(ynm1) + Labmda(ynm)) / 2;
            double lnm_12 = (Labmda(ynm_1) + Labmda(ynm)) / 2;
            double dynm_12 = (ynm_1 - ynm) / (hz * hz);
            double dynm12 = (ynm - ynm1) / (hz * hz);
            var t4 = (2 * ynm / tau) + (lnm12 * (ynm_1 - 2 * ynm + ynm1) / (hz * hz)) + F(x, z);
            return t4;
            // return 2 * ynm / tau  - lnm_12 * dynm_12 + lnm12 * dynm12 + F(x, z);
        }

        private double Am()
        {
            double lnm_12 = Labmda(1);
            return lnm_12 / (hz * hz);
        }

        private double Bm()
        {
            return Am() + Cm() + 2 / tau;
        }

        private double Cm()
        {
            double lnm12 = Labmda(1);
            return lnm12 / (hz * hz);
        }

        private double Dm(double yn_1m, double ynm, double yn1m, double x, double z)
        {
            double ln12m = (Labmda(yn1m) + Labmda(ynm)) / 2;
            double ln_12m = (Labmda(yn_1m) + Labmda(ynm)) / 2;
            double dyn_12m = (yn_1m - ynm) / (hx * hx);
            double dyn12m = (ynm - yn1m) / (hx * hx);
            return 2 * ynm / tau + ln12m * (yn_1m - 2 * ynm + yn1m) / (hx * hx) + F(x, z);
            //return 2 * ynm / tau - ln_12m * dyn_12m + ln12m * dyn12m + F(x, z);
        }

        public double[] RunMethodX(double[] a, double[] b, double[] c, double[] d)
        {
            double[] y = new double[b.Length];
            double[] xi = new double[b.Length - 1];
            double[] eta = new double[b.Length - 1];
            // xi[0] = -k0 / m0;// 0
            // eta[0] = p0 / m0;// u0
            xi[0] = 1;
            eta[0] = S * hx / l;
            for (int i = 1; i < xi.Length; i++)
            {
                double tmp = b[i] - a[i] * xi[i - 1];
                xi[i] = c[i] / tmp;
                eta[i] = (a[i] * eta[i - 1] + d[i]) / tmp;
            }
            //y[^1] = (pn - mn * eta[^1]) / (kn + mn * xi[^1]); //u0
            y[^1] = T0;//u0
            for (int i = y.Length - 1; i > 0; i--)
            {
                y[i - 1] = xi[i - 1] * y[i] + eta[i - 1];
            }
            return y;
        }

        public double[] RunMethodZ(double[] a, double[] b, double[] c, double[] d)
        {
            double[] y = new double[b.Length];
            double[] xi = new double[b.Length - 1];
            double[] eta = new double[b.Length - 1];
            //xi[0] = -k0 / m0;// 0
            //eta[0] = p0 / m0;// u0
            xi[0] = 1;
            eta[0] = 0;
            for (int i = 1; i < xi.Length; i++)
            {
                double tmp = b[i] - a[i] * xi[i - 1];
                xi[i] = c[i] / tmp;
                eta[i] = (a[i] * eta[i - 1] + d[i]) / tmp;
            }
            //y[^1] = (pn - mn * eta[^1]) / (kn + mn * xi[^1]); //u0
            //y[^1] = T0; //u0
            y[^1] = eta[^1] / (1 - xi[^1]);
            for (int i = y.Length - 1; i > 0; i--)
            {
                y[i - 1] = xi[i - 1] * y[i] + eta[i - 1];
            }
            return y;
        }

        public double[] RunMethod(double[] a, double[] b, double[] c, double[] d)
        {
            double[] y = new double[b.Length];
            double[] xi = new double[b.Length - 1];
            double[] eta = new double[b.Length - 1];
            // xi[0] = -k0 / m0;// 0
            // eta[0] = p0 / m0;// u0
            xi[0] = 0;
            eta[0] = T0;
            for (int i = 1; i < xi.Length; i++)
            {
                double tmp = b[i] - a[i] * xi[i - 1];
                xi[i] = c[i] / tmp;
                eta[i] = (a[i] * eta[i - 1] + d[i]) / tmp;
            }
            //y[^1] = (pn - mn * eta[^1]) / (kn + mn * xi[^1]); //u0
            y[^1] = T0; //u0
            for (int i = y.Length - 1; i > 0; i--)
            {
                y[i - 1] = xi[i - 1] * y[i] + eta[i - 1];
            }
            return y;
        }

        public double[] NextLayerX(double[] ym, double[] ym_1, double[] ym1, double[] x, double[] z,
            Double2Func A, Double2Func B, Double2Func C, Double5Func D)

        {
            int n = ym.Length;
            // double[] y = new double[n];
            // ym.CopyTo(y, 0);
            // double[] erry = new double[n];
            // double maxerr;
           //  int its = 0;

            double[] a = new double[n], b = new double[n], c = new double[n], d = new double[n];
            for (int i = 1; i < n - 1; i++)
            {
                a[i] = A();
                b[i] = B();
                c[i] = C();
                d[i] = D(ym_1[i], ym[i], ym1[i], x[i], z[i]);
            }
            double[] newy = RunMethod(a, b, c, d);
            return newy;
        }

        public double[] NextLayerZ(double[] ym, double[] ym_1, double[] ym1, double[] x, double[] z,
            Double2Func A, Double2Func B, Double2Func C, Double5Func D)

        {
            int n = ym.Length;
            double[] y = new double[n];
            ym.CopyTo(y, 0);
            double[] erry = new double[n];
            double maxerr;
            int its = 0;

            double[] a = new double[n], b = new double[n], c = new double[n], d = new double[n];
            for (int i = 1; i < n - 1; i++)
            {
                a[i] = A();
                b[i] = B();
                c[i] = C();
                d[i] = D(ym_1[i], ym[i], ym1[i], x[i], z[i]);
            }
            double[] newy = RunMethod(a, b, c, d);
            return newy;
        }

        public double[] NextLayer(double[] ym, double[] ym_1, double[] ym1, double[] x, double[] z,
            Double2Func A, Double2Func B, Double2Func C, Double5Func D)

        {
            int n = ym.Length;
            double[] y = new double[n];
            ym.CopyTo(y, 0);
            double[] erry = new double[n];
            double maxerr;
            int its = 0;

            double[] a = new double[n], b = new double[n], c = new double[n], d = new double[n];
            for (int i = 1; i < n - 1; i++)
            {
                a[i] = A();
                b[i] = B();
                c[i] = C();
                d[i] = D(ym_1[i], ym[i], ym1[i], x[i], z[i]);
            }
            double[] newy = RunMethod(a, b, c, d);
            return newy;
        }

        private static double[][] Transpose(double[][] m)
        {
            double[][] mt = new double[m[0].Length][];
            for (int i = 0; i < mt.Length; i++)
                mt[i] = new double[m.Length];
            for (int i = 0; i < m.Length; i++)
                for (int j = 0; j < m[i].Length; j++)
                    mt[j][i] = m[i][j];
            return mt;
        }

        private static double[] Double2Arr(double a, int size)
        {
            double[] res = new double[size];
            for (int i = 0; i < size; i++) res[i] = a;
            return res;
        }

        public double[][] NextTime(double[] x, double[] z, double[][] y_)
        {
            int n = x.Length, m = z.Length;
            double[][] midy = new double[m][];
            midy[0] = new double[n];

            for (int i = 1; i < m - 1; i++)
                midy[i] = NextLayerX(y_[i], y_[i - 1], y_[i + 1], x, Double2Arr(z[i], m), An, Bn, Cn, Dn);
            for (int i = 0; i < n; i++)
                midy[0][i] = midy[1][i];
            midy[^1] = new double[n];
            for (int i = 0; i < n; i++)
                midy[m - 1][i] = midy[m - 2][i];
            midy = Transpose(midy);
            double[][] newy = new double[n][];
            newy[0] = new double[m];
            for (int i = 1; i < n - 1; i++)
                newy[i] = NextLayerZ(midy[i], midy[i - 1], midy[i + 1], Double2Arr(x[i], n), z, Am, Bm, Cm, Dm);

            for (int i = 0; i < m; i++)
                newy[0][i] = newy[1][i];

            newy[^1] = new double[m];
            for (int i = 0; i < m; i++)
                newy[n - 1][i] = newy[n - 2][i];
            return Transpose(newy);
        }
    }
}
