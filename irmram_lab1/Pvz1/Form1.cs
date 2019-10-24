using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace Pvz1
{
    public partial class Form1 : Form
    {
        List<Timer> Timerlist = new List<Timer>();

        public Form1()
        {
            InitializeComponent();
            Initialize();
            FillIntervals();
            comboBox1.SelectedIndex = 0;

        }

        // ---------------------------------------------- STYGU METODAS ----------------------------------------------

        double x1, x2, xtemp, k; // izoliacijos intervalo pradžia ir galas, vidurio taškas, koeficientas
        int N = 100; // maksimalus iteracijų skaičius
        int iii; // iteracijos numeris
        private double stepSize = 0.1;
        double nearest = 0;

        Series Fx, X1X2, XMid, Gx, XY; // naudojama atvaizduoti f-jai, šaknų rėžiams ir vidiniams taškams

        /// <summary>
        /// Sprendžiama lygtis F(x) = 0
        /// </summary>
        /// <param name="x">funkcijos argumentas</param>
        /// <returns></returns>
        private double F(double x)
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    return (double)(-0.7 * Math.Pow(x, 4) + 4.16 * Math.Pow(x, 3) + 1.19 * Math.Pow(x, 2) - 33.4 * x + 31.51);
                case 1:
                    return (double)((double)(Math.Pow(Math.E, Math.Sin(x))) - (double)(x/10));
                case 2:
                    return (double)(Math.Pow(x, 3) - (3 * Math.Pow(x, 2)) + (double)(1.5 / Math.PI));
            }
            return 0;
        }

        private double G(double x)
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    return (double)(x - ((-0.7 * Math.Pow(x, 4) + 4.16 * Math.Pow(x, 3) + 1.19 * Math.Pow(x, 2) - 33.4 * x + 31.51) / (-2.8 * Math.Pow(x, 3) + 12.48 * Math.Pow(x, 2) + 2.38 * x - 33.4)));
                case 1:
                    return (double)(x - ((Math.Pow(Math.E, Math.Sin(x)) - x / 10) / (Math.Pow(Math.E, Math.Sin(x)) * Math.Cos(x) - 0.1)));
                case 2:
                    return (double)(x- ((Math.Pow(x, 3) - (3 * Math.Pow(x, 2)) + (double)(1.5 / Math.PI)) / (3 * Math.Pow(x, 2) - 6 * x)));
            }
            return 0;
        }

        //private double dG(double x)
        //{
        //    switch (comboBox1.SelectedIndex)
        //    {
        //        case 0:
        //            return (double)(1 - 
        //                ((1.96 * Math.Pow(x, 6)
        //                - 17.472 * Math.Pow(x, 5)
        //                + 50.2508 * Math.Pow(x, 4)
        //                - 73.7184 * Math.Pow(x, 3)
        //                + 267.5162 * Math.Pow(x, 2)
        //                - 865.9816 * x + 1040.5662) / 
        //                Math.Pow(-2.8 * Math.Pow(x, 3) + 12.48 * Math.Pow(x, 2) + 2.38 * x - 33.4, 2)
        //                ));
        //        case 1:
        //            return (double)(10 * Math.Pow(Math.E, Math.Sin(x)) * Math.Cos(x));
        //        case 2:
        //            return (double)(3 * Math.Pow(x, 2) - (6 * x) + 1);
        //    }
        //    return 0;
        //}

        private Queue<Tuple<double, double>> int1 = new Queue<Tuple<double, double>>();
        private Queue<Tuple<double, double>> int2 = new Queue<Tuple<double, double>>();
        private Queue<Tuple<double, double>> int3 = new Queue<Tuple<double, double>>();

        private Queue<Tuple<double, double>> intervals
        {
            get
            {
                switch (comboBox1.SelectedIndex)
                {
                    case 0: return int1;
                    case 1: return int2;
                    case 2: return int3;
                    default: return null;
                }
            }
        }

        private void FillIntervals()
        {
            int1.Clear();
            int2.Clear();
            int3.Clear();

            comboBox1.SelectedIndex = 0;
            x1 = -3;
            x2 = 5;
            var xTemp = x1;
            for (x1 += stepSize; x1 < x2; x1 += stepSize)
            {
                if (Math.Sign(F(x1)) != Math.Sign(F(xTemp)))
                {
                    Console.WriteLine("{0}  {1}",xTemp, x1);
                    int1.Enqueue(new Tuple<double, double>(xTemp, x1));
                }

                xTemp = x1;
            }

            comboBox1.SelectedIndex = 1;
            x1 = 1;
            x2 = 15;
            xTemp = x1;
            for (x1 += stepSize; x1 < x2; x1 += stepSize)
            {
                if (Math.Sign(F(x1)) != Math.Sign(F(xTemp)))
                {
                    Console.WriteLine("{0}  {1}", xTemp, x1);
                    int2.Enqueue(new Tuple<double, double>(xTemp, x1));
                }

                xTemp = x1;
            }

            comboBox1.SelectedIndex = 2;
            x1 = -1;
            x2 = 4;
            xTemp = x1;
            for (x1 += stepSize; x1 < x2; x1 += stepSize)
            {
                if (Math.Sign(F(x1)) != Math.Sign(F(xTemp)))
                {
                    Console.WriteLine("{0}  {1}", xTemp, x1);
                    int3.Enqueue(new Tuple<double, double>(xTemp, x1));
                }

                xTemp = x1;
            }
        }

        // Mygtukas "Pusiaukirtos metodas" - ieškoma šaknies, ir vizualizuojamas paieškos procesas
        private void button3_Click(object sender, EventArgs e)
        {
            ClearForm(); // išvalomi programos duomenys
            prepareForm();
            x1 = -3; // izoliacijos intervalo pradžia
            x2 = -2; // izoliacijos intervalo galas
            iii = 0; // iteraciju skaičius
            richTextBox1.AppendText("Iteracija         x            F(x)        x1          x2          F(x1)         F(x2)       \n");
            // Nubraižoma f-ja, kuriai ieskome saknies
            Fx = chart1.Series.Add("F(x)");
            Fx.ChartType = SeriesChartType.Line;
            double x = -8;
            for (int i = 0; i < 300; i++)
            {
                Fx.Points.AddXY(x, F(x));
                x = x + 0.1;
            }
            Fx.BorderWidth = 3;
            
            X1X2 = chart1.Series.Add("X1X2");
            X1X2.MarkerStyle = MarkerStyle.Circle;
            X1X2.MarkerSize = 8;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;


            XMid = chart1.Series.Add("XMid");
            XMid.MarkerStyle = MarkerStyle.Circle;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;
            XMid.MarkerSize = 8;

            timer2.Enabled = true;
            timer2.Interval = 500; // timer2 intervalas milisekundemis
            timer2.Start();           
        }
        

        /// <summary>
        /// timer2 iteracijoje atliekami veiksmai
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void timer2_Tick(object sender, EventArgs e)
        {
            k = (float)(Math.Abs(F(x1) / F(x2)));
            xtemp = (x1 + k * x2) / (1 + k);
            //xtemp = (x1 + x2) / 2; // apskaiciuojamas vidurinis taskas


            if (Math.Abs(F(xtemp)) > 1e-4 & iii <= N)
            // tikrinama salyga, ar funkcijos absoliuti reiksme daugiau uz nustatyta (norima) 
            // tiksluma ir nevirsytas maksimalus iteraciju skaicius
            {
                X1X2.Points.Clear();
                XMid.Points.Clear();

                X1X2.Points.AddXY(x1, 0);
                X1X2.Points.AddXY(x2, 0);
                XMid.Points.AddXY(xtemp, 0);

                richTextBox1.AppendText(String.Format(" {0,6:d}   {1,12:f7}  {2,12:f7} {3,12:f7} {4,12:f7} {5,12:f7} {6,12:f7}\n",
                    iii, xtemp, F(xtemp), x1, x2, F(x1), F(x2)));
                if (Math.Sign((double)F(x1)) != Math.Sign((double)F(xtemp)))
                {
                    x2 = xtemp;
                }
                else
                {
                    x1 = xtemp;
                }
                iii = iii + 1;

            }
            else
            // skaiciavimai stabdomi
            {
                richTextBox1.AppendText("Skaičiavimai baigti");
                timer2.Stop();
            }
        }


        // ---------------------------------------------- PARAMETRINĖS FUNKCIJOS ----------------------------------------------

        List<PointF> data = new List<PointF>(); 
        Series S1;

        /// <summary>
        /// Parametrinis interpoliavimas
        /// </summary>
        private void button5_Click(object sender, EventArgs e)
        {
            ClearForm(); // išvalomi programos duomenys
            prepareForm();
            data.Clear();
            // apskaičiuojamos funkcijos reikšmės
            for (int i = 0; i < 400; i++)
            {
                float x = i / 50f * (float)(Math.Sin(2*i / 10f));
                float y = i / 50f * (float)(Math.Sin(i / 10f));
                data.Add(new PointF(x, y));
            }
            S1 = chart1.Series.Add("S1");
            S1.BorderWidth = 9;
            S1.ChartType = SeriesChartType.Line;
          
            timer3.Enabled = true;
            timer3.Interval = 15;
            timer3.Start();
        }

        private void timer3_Tick(object sender, EventArgs e)
        {
            Series S1 = chart1.Series[0];
            int pointsSoFar = S1.Points.Count;
            if (pointsSoFar < data.Count)
            {
                S1.Points.AddXY(data[pointsSoFar].X, data[pointsSoFar].Y);
            }
            else
            {
                timer1.Stop();
            }
        }


        // ---------------------------------------------- TIESINĖ ALGEBRA ----------------------------------------------

        /// <summary>
        /// Tiesine algebra (naudojama MathNet)
        /// </summary>
        private void button2_Click(object sender, EventArgs e)
        {
            ClearForm();

            double[,] x = { { 1, 2, 3 }, { 3, 4, 5 }, { 6, 5, 8 } };
            // iš masyvo sugeneruoja matricą, is matricos išskiria eilutę - suformuoja vektorių
            Matrix<double> m = Matrix<double>.Build.DenseOfArray(x);
            Vector<double> v = m.Row(1);
            richTextBox1.AppendText("\nMatrica m:\n");
            richTextBox1.AppendText(m.ToString());

            richTextBox1.AppendText("\nVektorius v:\n");
            richTextBox1.AppendText(v.ToString());

            richTextBox1.AppendText("\ntranspose(m):\n");
            richTextBox1.AppendText(m.Transpose().ToString());

            Matrix<double> vm = v.ToRowMatrix();
            richTextBox1.AppendText("\nvm = v' - toRowMatrix()\n");
            richTextBox1.AppendText(vm.ToString());

            Vector<double> v1 = m * v;
            richTextBox1.AppendText("\nv1 = m * v\n");
            richTextBox1.AppendText(v1.ToString());
            richTextBox1.AppendText("\nmin(v1)\n");
            richTextBox1.AppendText(v1.Min().ToString());

            Matrix<double> m1 = m.Inverse();
            richTextBox1.AppendText("\ninverse(m)\n");
            richTextBox1.AppendText(m1.ToString());

            richTextBox1.AppendText("\ndet(m)\n");
            richTextBox1.AppendText(m.Determinant().ToString());

            // you must add reference to assembly system.Numerics
            Evd<double> eigenv = m.Evd();
            richTextBox1.AppendText("\neigenvalues(m)\n");
            richTextBox1.AppendText(eigenv.EigenValues.ToString());
            
            LU<double> LUanswer = m.LU();
            richTextBox1.AppendText("\nMatricos M LU skaida\n");
            richTextBox1.AppendText("\nMatrica L:\n");
            richTextBox1.AppendText(LUanswer.L.ToString());
            richTextBox1.AppendText("\nMatrica U:\n");
            richTextBox1.AppendText(LUanswer.U.ToString());
            
            QR<double> QRanswer = m.QR();
            richTextBox1.AppendText("\nMatricos M QR skaida\n");
            richTextBox1.AppendText("\nMatrica Q:\n");
            richTextBox1.AppendText(QRanswer.Q.ToString());
            richTextBox1.AppendText("\nMatrica R:\n");
            richTextBox1.AppendText(QRanswer.R.ToString());

            Vector<double> v3 = m.Solve(v);
            richTextBox1.AppendText("\nm*v3 = v sprendziama QR metodu\n");
            richTextBox1.AppendText(v3.ToString());
            richTextBox1.AppendText("Patikrinimas\n");
            richTextBox1.AppendText((m * v3 - v).ToString());
            
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        // ---------------------------------------------- PAPRASTUJU ITERACIJU METODAS ------------------------------

        private void button7_Click(object sender, EventArgs e)
        {
            ClearForm();

            double x = -8;
            int i = 0;
            prepareForm();
            switch (comboBox1.SelectedIndex)
            {
                case 2:
                    i = 1;
                    break;
            }

            richTextBox1.AppendText("Iteracija         x            F(x)        x1          x2          F(x1)         F(x2)       \n");

            Fx = chart1.Series.Add("F(x)");
            Fx.ChartType = SeriesChartType.Line;
            Gx = chart1.Series.Add("G(x)");
            Gx.ChartType = SeriesChartType.Line;
            XY = chart1.Series.Add("y = x");
            XY.ChartType = SeriesChartType.Line;
            for (; i < 300; i++)
            {
                Fx.Points.AddXY(x, F(x));
                Gx.Points.AddXY(x, G(x));
                XY.Points.AddXY(x, x); x = x + (2 * Math.PI) / 50;
            }
            Fx.BorderWidth = 3;
            Gx.BorderWidth = 3;
            XY.BorderWidth = 3;

            X1X2 = chart1.Series.Add("X1X2");
            X1X2.MarkerStyle = MarkerStyle.Circle;
            X1X2.MarkerSize = 8;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;


            XMid = chart1.Series.Add("XMid");
            XMid.MarkerStyle = MarkerStyle.Circle;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;
            XMid.MarkerSize = 8;

            var thing = intervals.Dequeue();
            x1 = thing.Item1;

            timer6.Enabled = true;
            timer6.Interval = 500;
            timer6.Start();
        }

        private void timer6_Tick(object sender, EventArgs e)
        {
            xtemp = G(x1);
            x1 = xtemp;

            if (Math.Abs(F(x1)) > 1e-6 & iii <= N)
            {
                X1X2.Points.Clear();
                XMid.Points.Clear();

                XMid.Points.AddXY(x1, 0);

                richTextBox1.AppendText($" {iii,6:d}   {xtemp,12:f7}  {F(xtemp),12:f7} {x1,12:f7} {x2,12:f7} {F(x1),12:f7} {F(x2),12:f7}\n");

                iii = iii + 1;
            }
            else
            {
                richTextBox1.AppendText($" {iii,6:d}   {xtemp,12:f7}  {F(xtemp),12:f7} {x1,12:f7} {x2,12:f7} {F(x1),12:f7} {F(x2),12:f7}\n");
                richTextBox1.AppendText("Skaičiavimai baigti");
                iii = 0;
                if (intervals.Any())
                {
                    var thing = intervals.Dequeue();
                    x1 = thing.Item1;
                }
                else
                {
                    timer6.Stop();
                }
            }
        }

        // ---------------------------------------------- INTERVALU SKENAVIMO METODAS -------------------------------

        private void button8_Click(object sender, EventArgs e)
        {
            ClearForm();

            double x = -8;
            int i = 0;
            if (comboBox1.SelectedIndex == 2)
            {
                //x = 0.01;
                i = 1;
            }
            prepareForm();

            richTextBox1.AppendText("Iteracija         x            F(x)        x1          x2          F(x1)         F(x2)       \n");

            Fx = chart1.Series.Add("F(x)");
            Fx.ChartType = SeriesChartType.Line;
            for (; i < 300; i++)
            {
                Fx.Points.AddXY(x, F(x)); x = x + (2 * Math.PI) / 50;
            }
            Fx.BorderWidth = 3;

            X1X2 = chart1.Series.Add("X1X2");
            X1X2.MarkerStyle = MarkerStyle.Circle;
            X1X2.MarkerSize = 8;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;


            XMid = chart1.Series.Add("XMid");
            XMid.MarkerStyle = MarkerStyle.Circle;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;
            XMid.MarkerSize = 8;

            var thing = intervals.Dequeue();
            x1 = thing.Item1;
            x2 = thing.Item2;

            timer5.Enabled = true;
            timer5.Interval = 50;
            timer5.Start();
        }

        private void timer5_Tick(object sender, EventArgs e)
        {
            xtemp = x1 + stepSize;

            if (Math.Abs(F(xtemp)) > 1e-6 & iii <= N)
            {
                X1X2.Points.Clear();

                X1X2.Points.AddXY(x1, 0);
                X1X2.Points.AddXY(xtemp, 0);

                richTextBox1.AppendText($" {iii,6:d}   {xtemp,12:f7}  {F(xtemp),12:f7} {x1,12:f7} {x2,12:f7} {F(x1),12:f7} {F(x2),12:f7}\n");

                if (Math.Sign(F(x1)) != Math.Sign(F(xtemp)))
                {
                    stepSize /= 10;
                }
                else
                {
                    x1 += stepSize;
                }


                iii = iii + 1;
            }
            else
            {
                richTextBox1.AppendText($" {iii,6:d}   {xtemp,12:f7}  {F(xtemp),12:f7} {x1,12:f7} {x2,12:f7} {F(x1),12:f7} {F(x2),12:f7}\n");
                richTextBox1.AppendText("Skaičiavimai baigti");
                iii = 0;
                if (intervals.Any())
                {
                    var thing = intervals.Dequeue();
                    x1 = thing.Item1;
                    x2 = thing.Item2;
                    stepSize = 0.1;
                }
                else
                {
                    timer5.Stop();
                }
            }
        }

        // ---------------------------------------------- STYGU METODAS ---------------------------------------------

        private void Button6_Click(object sender, EventArgs e)
        {
            ClearForm(); // išvalomi programos duomenys
            prepareForm();
            x1 = -3; // izoliacijos intervalo pradžia
            x2 = -2; // izoliacijos intervalo galas
            iii = 0; // iteraciju skaičius
            richTextBox1.AppendText("Iteracija         x            F(x)        x1          x2          F(x1)         F(x2)       \n");
            // Nubraižoma f-ja, kuriai ieskome saknies
            Fx = chart1.Series.Add("F(x)");
            Fx.ChartType = SeriesChartType.Line;
            double x = -8;
            for (int i = 0; i < 300; i++)
            {
                Fx.Points.AddXY(x, F(x));
                x = x + 0.1;
            }
            Fx.BorderWidth = 3;

            X1X2 = chart1.Series.Add("X1X2");
            X1X2.MarkerStyle = MarkerStyle.Circle;
            X1X2.MarkerSize = 8;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;


            XMid = chart1.Series.Add("XMid");
            XMid.MarkerStyle = MarkerStyle.Circle;
            X1X2.ChartType = SeriesChartType.Point;
            X1X2.ChartType = SeriesChartType.Line;
            XMid.MarkerSize = 8;

            var thing = intervals.Dequeue();
            x1 = thing.Item1;
            x2 = thing.Item2;

            timer4.Enabled = true;
            timer4.Interval = 500; // timer2 intervalas milisekundemis
            timer4.Start();
        }

        private void Timer4_Tick(object sender, EventArgs e)
        {
            k = (float)(Math.Abs(F(x1) / F(x2)));
            xtemp = (x1 + k * x2) / (1 + k);

            if (Math.Abs(F(xtemp)) > 1e-6 & iii <= N)
            // tikrinama salyga, ar funkcijos absoliuti reiksme daugiau uz nustatyta (norima) 
            // tiksluma ir nevirsytas maksimalus iteraciju skaicius
            {
                X1X2.Points.Clear();
                XMid.Points.Clear();

                X1X2.Points.AddXY(x1, F(x1));
                X1X2.Points.AddXY(x2, F(x2));
                XMid.Points.AddXY(xtemp, F(xtemp));

                richTextBox1.AppendText(String.Format(" {0,6:d}   {1,12:f7}  {2,12:f7} {3,12:f7} {4,12:f7} {5,12:f7} {6,12:f7}\n",
                    ++iii, xtemp, F(xtemp), x1, x2, F(x1), F(x2)));
                if (Math.Sign((double)F(x1)) != Math.Sign((double)F(xtemp)))
                {
                    x2 = xtemp;
                }
                else
                {
                    x1 = xtemp;
                }

            }
            else
            {
                richTextBox1.AppendText($" {++iii,6:d}   {xtemp,12:f7}  {F(xtemp),12:f7} {x1,12:f7} {x2,12:f7} {F(x1),12:f7} {F(x2),12:f7}\n");
                richTextBox1.AppendText("Skaičiavimai baigti\n");
                iii = 0;
                if (intervals.Any())
                {
                    var thing = intervals.Dequeue();
                    x1 = thing.Item1;
                    x2 = thing.Item2;
                }
                else
                {
                    timer4.Stop();
                }
            }
        }


        // ---------------------------------------------- KITI METODAI ----------------------------------------------

        /// <summary>
        /// Uždaroma programa
        /// </summary>
        private void button1_Click(object sender, EventArgs e)
        {
            Close();
        }
        
        /// <summary>
        /// Išvalomas grafikas ir consolė
        /// </summary>
        private void button4_Click(object sender, EventArgs e)
        {
            ClearForm();
        }
        

        public void ClearForm()
        {
            richTextBox1.Clear(); // isvalomas richTextBox1
            // sustabdomi timeriai jei tokiu yra
            foreach (var timer in Timerlist)
            {
                timer.Stop();
            }

            // isvalomos visos nubreztos kreives
            var s = comboBox1.SelectedIndex;
            FillIntervals();
            comboBox1.SelectedIndex = s;

            stepSize = 0.1;
            richTextBox1.Clear();
            chart1.Series.Clear();
        }

        /// <summary>
        /// Paruosia grafiko forma atitinkamai funkcijai
        /// </summary>
        private void prepareForm()
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    PreparareForm(-5, 5, -5, 5);
                    break;
                case 1:
                    PreparareForm(1, 15, -3, 3);
                    break;
                case 2:
                    PreparareForm(-4, 4, -4, 4);
                    break;
            }
        }
    }
}
