using System.Runtime.CompilerServices;
using lab_05;

namespace lab_05
{
    internal class Program
    {
        static void Main(string[] args)
        {
            double xmax = 10, zmax = 10;
            int n = 5, m = 5;
            double[] x = new double[n], z = new double[m];
            for (int i = 0; i < n; i++)
                x[i] = i * xmax / n;
            for (int i = 0; i < m; i++)
                z[i] = i * zmax / m;
            var mod = new Model(
                a1: 0.0134, 
                b1: 1, 
                c1: 4.35e-4, 
                m1: 1, 
                f0: 300, 
                beta: 0.5,
                x0: 5, 
                z0: 5, 
                hx: xmax / n, hz: zmax / m, 
                F0: 0, 
                T0: 300, 
                alpha2: 0.05,
                alpha3: 0.05, 
                alpha4: 0.05);
            double[][] y = new double[m][];
            for (int i = 0; i < n; i++)
            {
                y[i] = new double[n];
                for (int j = 0; j < m; j++)
                    y[i][j] = mod.T0;
            }
            double maxerr;
            double[][] materr = new double[m][];
            for (int i = 0; i < n; i++)
                materr[i] = new double[n];
            double[] maxerra = new double[m];
            int its = 0;
            do
            {
                double[][] newy = mod.NextTime(x, z, y);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                        materr[i][j] = Math.Abs((newy[i][j] - y[i][j]) / newy[i][j]);
                    maxerra[i] = materr[i].Max();
                }
                maxerr = maxerra.Max();
                y = newy;
                its++;
                Console.WriteLine(its);
            } while (maxerr > 1e-5);

            using (StreamWriter sw = new StreamWriter("data.txt"))
            {
                sw.WriteLine($"{n} {m}");
                sw.WriteLine(string.Join(" ", x));
                sw.WriteLine(string.Join(" ", z));
                for (int i = 0; i < m; i++)
                    sw.WriteLine(string.Join(" ", y[i]));
            }

        }
    }
}