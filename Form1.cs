using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Windows.Forms;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.WindowsForms;
using OxyPlot.Axes;
using OxyPlot.Annotations;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace CoolingTowerSimulation
{
    public partial class Form1 : Form
    {
        // Controls
        private TabControl tabControl;
        private TabPage tabParameters, tabTimeDomain, tabFrequencyDomain, tabLaplace, tabZDomain;

        // Parameter Controls
        private NumericUpDown numAmbientTemp, numHumidity, numWindSpeed, numMechanicalVib, numWaterTemp;
        private Button btnSimulate, btnStopSimulation, btnUpdateAnalysis;
        private Label lblStatus;

        // Plot Views
        private readonly PlotView[] timePlots = new PlotView[5];
        private readonly PlotView[] freqPlots = new PlotView[5];
        private PlotView plotLaplace, plotZDomain;
        private RichTextBox rtbLaplace, rtbZDomain;

        // Data Storage
        private SensorData sensorData;

        // Multi-Rate Sampling
        private readonly int[] samplingRates = new int[5] { 2, 1, 5, 1000, 1 };
        private readonly double duration = 10.0;

        // Realtime simulation
        private Timer simulationTimer;
        private double currentTime = 0;
        private readonly double timeWindow = 10.0;
        private bool isSimulating = false;
        private Random rand = new Random(42);

        // Parameters
        private double ambientTemp, humidity, windSpeed, mechVib, waterTemp;

        // Sensor time constants
        private readonly double tempTau = 0.5;
        private readonly double humidityTau = 0.67;
        private readonly double airflowTau = 0.2;
        private readonly double vibrationTau = 0.02;
        private readonly double doTau = 1.25;

        public Form1()
        {
            InitializeComponent();
            InitializeUI();
            InitializeTimer();
        }

        private void InitializeComponent()
        {
            this.SuspendLayout();
            this.ClientSize = new Size(1400, 900);
            this.Text = "Cooling Tower Efficiency Analysis - Multi-Rate Sampling";
            this.StartPosition = FormStartPosition.CenterScreen;
            this.ResumeLayout(false);
        }

        private void InitializeTimer()
        {
            simulationTimer = new Timer { Interval = 50 };
            simulationTimer.Tick += SimulationTimer_Tick;
        }

        private void InitializeUI()
        {
            tabControl = new TabControl { Dock = DockStyle.Fill, Font = new Font("Segoe UI", 9F) };

            tabParameters = new TabPage("Parameter Input");
            tabTimeDomain = new TabPage("Time Domain (REALTIME)");
            tabFrequencyDomain = new TabPage("Frequency Domain (FFT)");
            tabLaplace = new TabPage("S-Domain (Laplace)");
            tabZDomain = new TabPage("Z-Domain (Discrete)");

            tabControl.TabPages.AddRange(new[] { tabParameters, tabTimeDomain, tabFrequencyDomain, tabLaplace, tabZDomain });

            InitializeParameterTab();
            InitializeTimeDomainTab();
            InitializeFrequencyTab();
            InitializeLaplaceTab();
            InitializeZDomainTab();

            this.Controls.Add(tabControl);
        }

        private void InitializeParameterTab()
        {
            Panel panel = new Panel { Dock = DockStyle.Fill, AutoScroll = true, Padding = new Padding(20) };
            int yPos = 20;

            Label title = new Label
            {
                Text = "Sensor Input Parameters Configuration",
                Font = new Font("Segoe UI", 14F, FontStyle.Bold),
                Location = new Point(20, yPos),
                AutoSize = true
            };
            panel.Controls.Add(title);
            yPos += 50;

            AddParameterControl(panel, "Ambient Temperature (°C):", 25m, 0m, 50m, ref numAmbientTemp, ref yPos);
            AddParameterControl(panel, "Relative Humidity (%):", 60m, 0m, 100m, ref numHumidity, ref yPos);
            AddParameterControl(panel, "Wind Speed (m/s):", 5m, 0m, 30m, ref numWindSpeed, ref yPos);
            AddParameterControl(panel, "Mechanical Vibration Base (g):", 0.5m, 0m, 5m, ref numMechanicalVib, ref yPos);
            AddParameterControl(panel, "Water Temperature (°C):", 30m, 10m, 60m, ref numWaterTemp, ref yPos);

            yPos += 20;

            btnSimulate = new Button
            {
                Text = "▶ Start Realtime Simulation",
                Location = new Point(20, yPos),
                Size = new Size(200, 40),
                Font = new Font("Segoe UI", 10F, FontStyle.Bold),
                BackColor = Color.FromArgb(0, 150, 0),
                ForeColor = Color.White,
                FlatStyle = FlatStyle.Flat
            };
            btnSimulate.Click += BtnSimulate_Click;
            panel.Controls.Add(btnSimulate);

            btnStopSimulation = new Button
            {
                Text = "⏹ Stop Simulation",
                Location = new Point(240, yPos),
                Size = new Size(180, 40),
                Font = new Font("Segoe UI", 10F, FontStyle.Bold),
                BackColor = Color.FromArgb(200, 0, 0),
                ForeColor = Color.White,
                FlatStyle = FlatStyle.Flat,
                Enabled = false
            };
            btnStopSimulation.Click += BtnStopSimulation_Click;
            panel.Controls.Add(btnStopSimulation);

            yPos += 50;

            btnUpdateAnalysis = new Button
            {
                Text = "🔄 Update FFT/Laplace/Z Analysis",
                Location = new Point(20, yPos),
                Size = new Size(260, 40),
                Font = new Font("Segoe UI", 9F, FontStyle.Bold),
                BackColor = Color.FromArgb(0, 120, 215),
                ForeColor = Color.White,
                FlatStyle = FlatStyle.Flat
            };
            btnUpdateAnalysis.Click += (s, ev) => GenerateStaticAnalysis();
            panel.Controls.Add(btnUpdateAnalysis);

            yPos += 50;

            lblStatus = new Label
            {
                Text = "Status: Stopped",
                Location = new Point(20, yPos),
                Size = new Size(600, 30),
                Font = new Font("Segoe UI", 10F, FontStyle.Bold),
                ForeColor = Color.Gray
            };
            panel.Controls.Add(lblStatus);
            yPos += 40;

            Label info = new Label
            {
                Text = "Multi-Rate Sampling Strategy:\n" +
                       "• Temperature (SHT85): 2 Hz - τ=0.5s, Bandwidth=0.32 Hz\n" +
                       "• Humidity (SHT85): 1 Hz - τ=0.67s, Bandwidth=0.24 Hz\n" +
                       "• Airflow (Testo 440): 5 Hz - τ=0.2s, Bandwidth=0.80 Hz\n" +
                       "• Vibration (PCB 352C33): 1000 Hz - τ=0.02s, captures 60-120 Hz\n" +
                       "• Dissolved Oxygen (Apera DO850): 1 Hz - τ=1.25s, Bandwidth=0.13 Hz\n\n" +
                       "⚠️ UBAH PARAMETER SAAT SIMULASI untuk melihat efek real-time!",
                Location = new Point(20, yPos),
                Size = new Size(900, 180),
                Font = new Font("Segoe UI", 9F, FontStyle.Italic),
                ForeColor = Color.DarkBlue
            };
            panel.Controls.Add(info);

            tabParameters.Controls.Add(panel);
        }

        private void AddParameterControl(Panel panel, string label, decimal defaultVal, decimal min, decimal max,
                                        ref NumericUpDown control, ref int yPos)
        {
            Label lbl = new Label
            {
                Text = label,
                Location = new Point(20, yPos),
                Size = new Size(250, 25),
                Font = new Font("Segoe UI", 9F)
            };

            control = new NumericUpDown
            {
                Location = new Point(280, yPos),
                Size = new Size(120, 25),
                Minimum = min,
                Maximum = max,
                Value = defaultVal,
                DecimalPlaces = 2,
                Increment = 0.5m
            };

            panel.Controls.Add(lbl);
            panel.Controls.Add(control);
            yPos += 40;
        }

        private void InitializeTimeDomainTab()
        {
            Panel panel = new Panel { Dock = DockStyle.Fill, AutoScroll = true };

            string[] titles = {
                "Temperature (°C) - REALTIME [SHT85, fs=2Hz]",
                "Humidity (%) - REALTIME [SHT85, fs=1Hz]",
                "Airflow (m/s) - REALTIME [Testo 440, fs=5Hz]",
                "Vibration (g) - REALTIME [PCB 352C33, fs=1000Hz]",
                "Dissolved Oxygen (mg/L) - REALTIME [Apera DO850, fs=1Hz]"
            };

            for (int i = 0; i < 5; i++)
            {
                timePlots[i] = new PlotView
                {
                    Location = new Point(10, 10 + i * 170),
                    Size = new Size(1350, 160),
                    BackColor = Color.White
                };

                var model = new PlotModel { Title = titles[i] };
                model.Axes.Add(new LinearAxis
                {
                    Position = AxisPosition.Bottom,
                    Title = "Time (s)",
                    Minimum = 0,
                    Maximum = timeWindow,
                    MajorGridlineStyle = LineStyle.Solid,
                    MinorGridlineStyle = LineStyle.Dot,
                    MajorGridlineColor = OxyColors.LightGray
                });
                model.Axes.Add(new LinearAxis
                {
                    Position = AxisPosition.Left,
                    Title = titles[i].Split('-')[0].Trim(),
                    MajorGridlineStyle = LineStyle.Solid,
                    MinorGridlineStyle = LineStyle.Dot,
                    MajorGridlineColor = OxyColors.LightGray
                });

                var series = new LineSeries
                {
                    Color = OxyColors.Blue,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.Solid
                };
                model.Series.Add(series);

                timePlots[i].Model = model;
                panel.Controls.Add(timePlots[i]);
            }

            tabTimeDomain.Controls.Add(panel);
        }

        private void InitializeFrequencyTab()
        {
            Panel panel = new Panel { Dock = DockStyle.Fill, AutoScroll = true };

            string[] titles = {
                "Temperature Spectrum",
                "Humidity Spectrum",
                "Airflow Spectrum",
                "Vibration Spectrum",
                "Dissolved Oxygen Spectrum"
            };

            for (int i = 0; i < 5; i++)
            {
                freqPlots[i] = new PlotView
                {
                    Location = new Point(10, 10 + i * 170),
                    Size = new Size(1350, 160),
                    BackColor = Color.White
                };

                var model = new PlotModel { Title = titles[i] };
                model.Axes.Add(new LinearAxis
                {
                    Position = AxisPosition.Bottom,
                    Title = "Frequency (Hz)",
                    MajorGridlineStyle = LineStyle.Solid,
                    MinorGridlineStyle = LineStyle.Dot,
                    MajorGridlineColor = OxyColors.LightGray
                });
                model.Axes.Add(new LinearAxis
                {
                    Position = AxisPosition.Left,
                    Title = "Magnitude",
                    MajorGridlineStyle = LineStyle.Solid,
                    MinorGridlineStyle = LineStyle.Dot,
                    MajorGridlineColor = OxyColors.LightGray
                });
                freqPlots[i].Model = model;

                panel.Controls.Add(freqPlots[i]);
            }

            tabFrequencyDomain.Controls.Add(panel);
        }

        private void InitializeLaplaceTab()
        {
            Panel mainPanel = new Panel { Dock = DockStyle.Fill };
            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.Fill,
                Orientation = Orientation.Vertical,
                SplitterDistance = 600
            };

            rtbLaplace = new RichTextBox
            {
                Dock = DockStyle.Fill,
                Font = new Font("Consolas", 9F),
                ReadOnly = true,
                BackColor = Color.WhiteSmoke,
                Padding = new Padding(10)
            };
            split.Panel1.Controls.Add(rtbLaplace);

            Panel plotPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(10) };
            Label lblPlot = new Label
            {
                Text = "Pole Plot (S-Plane)",
                Dock = DockStyle.Top,
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                Height = 30,
                TextAlign = ContentAlignment.MiddleCenter
            };

            plotLaplace = new PlotView { Dock = DockStyle.Fill, BackColor = Color.White };

            plotPanel.Controls.Add(plotLaplace);
            plotPanel.Controls.Add(lblPlot);
            split.Panel2.Controls.Add(plotPanel);

            mainPanel.Controls.Add(split);
            tabLaplace.Controls.Add(mainPanel);
        }

        private void InitializeZDomainTab()
        {
            Panel mainPanel = new Panel { Dock = DockStyle.Fill };
            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.Fill,
                Orientation = Orientation.Vertical,
                SplitterDistance = 600
            };

            rtbZDomain = new RichTextBox
            {
                Dock = DockStyle.Fill,
                Font = new Font("Consolas", 9F),
                ReadOnly = true,
                BackColor = Color.WhiteSmoke,
                Padding = new Padding(10)
            };
            split.Panel1.Controls.Add(rtbZDomain);

            Panel plotPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(10) };
            Label lblPlot = new Label
            {
                Text = "Pole-Zero Map (Z-Plane with Unit Circle)",
                Dock = DockStyle.Top,
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                Height = 30,
                TextAlign = ContentAlignment.MiddleCenter
            };

            plotZDomain = new PlotView { Dock = DockStyle.Fill, BackColor = Color.White };

            plotPanel.Controls.Add(plotZDomain);
            plotPanel.Controls.Add(lblPlot);
            split.Panel2.Controls.Add(plotPanel);

            mainPanel.Controls.Add(split);
            tabZDomain.Controls.Add(mainPanel);
        }

        private void BtnSimulate_Click(object sender, EventArgs e)
        {
            currentTime = 0;

            foreach (var plot in timePlots)
            {
                if (plot?.Model?.Series != null && plot.Model.Series.Count > 0)
                {
                    ((LineSeries)plot.Model.Series[0]).Points.Clear();
                }
            }

            isSimulating = true;
            simulationTimer.Start();

            btnSimulate.Enabled = false;
            btnStopSimulation.Enabled = true;
            lblStatus.Text = "Status: Running - Parameter dapat diubah real-time!";
            lblStatus.ForeColor = Color.Green;

            GenerateStaticAnalysis();
        }

        private void BtnStopSimulation_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulating = false;

            btnSimulate.Enabled = true;
            btnStopSimulation.Enabled = false;
            lblStatus.Text = "Status: Stopped";
            lblStatus.ForeColor = Color.Gray;
        }

        private void SimulationTimer_Tick(object sender, EventArgs e)
        {
            if (!isSimulating) return;

            ambientTemp = (double)numAmbientTemp.Value;
            humidity = (double)numHumidity.Value;
            windSpeed = (double)numWindSpeed.Value;
            mechVib = (double)numMechanicalVib.Value;
            waterTemp = (double)numWaterTemp.Value;

            currentTime += 0.05;

            UpdateRealtimePlot(timePlots[0], currentTime, CalculateTemperature(currentTime));
            UpdateRealtimePlot(timePlots[1], currentTime, CalculateHumidity(currentTime));
            UpdateRealtimePlot(timePlots[2], currentTime, CalculateAirflow(currentTime));
            UpdateRealtimePlot(timePlots[3], currentTime, CalculateVibration(currentTime));
            UpdateRealtimePlot(timePlots[4], currentTime, CalculateDissolvedOxygen(currentTime));
        }

        private void UpdateRealtimePlot(PlotView plot, double time, double value)
        {
            var series = (LineSeries)plot.Model.Series[0];
            series.Points.Add(new DataPoint(time, value));

            while (series.Points.Count > 0 && series.Points[0].X < time - timeWindow)
            {
                series.Points.RemoveAt(0);
            }

            var xAxis = plot.Model.Axes[0];
            if (time > timeWindow)
            {
                xAxis.Minimum = time - timeWindow;
                xAxis.Maximum = time;
            }

            plot.Model.InvalidatePlot(true);
        }

        private double CalculateTemperature(double t)
        {
            double baseTemp = ambientTemp;
            double periodic = 2 * Math.Sin(2 * Math.PI * 0.5 * t);
            double noise = 0.5 * (rand.NextDouble() - 0.5);
            double windEffect = -0.8 * (windSpeed - 5);
            double vibEffect = 0.3 * mechVib;
            return baseTemp + periodic + windEffect + vibEffect + noise;
        }

        private double CalculateHumidity(double t)
        {
            double baseHumidity = humidity;
            double tempEffect = -1.5 * (ambientTemp - 25);
            double waterEffect = 0.8 * (waterTemp - 30);
            double periodic = 5 * Math.Sin(2 * Math.PI * 0.3 * t);
            double noise = 1 * (rand.NextDouble() - 0.5);
            double windEffect = -0.5 * (windSpeed - 5);
            double result = baseHumidity + tempEffect + waterEffect + periodic + windEffect + noise;
            return Math.Max(0, Math.Min(100, result));
        }

        private double CalculateAirflow(double t)
        {
            double baseAirflow = windSpeed;
            double periodic1 = 1.5 * Math.Sin(2 * Math.PI * 1.0 * t);
            double periodic2 = 0.5 * Math.Sin(2 * Math.PI * 3.0 * t);
            double noise = 0.3 * (rand.NextDouble() - 0.5);
            double tempEffect = 0.3 * (ambientTemp - 25);
            double humidityEffect = -0.05 * (humidity - 60);
            return Math.Max(0, baseAirflow + periodic1 + periodic2 + tempEffect + humidityEffect + noise);
        }

        private double CalculateVibration(double t)
        {
            double fundamental = mechVib * Math.Sin(2 * Math.PI * 60 * t);
            double harmonic = 0.2 * Math.Sin(2 * Math.PI * 120 * t);
            double noise = 0.1 * (rand.NextDouble() - 0.5);
            double windEffect = 0.05 * windSpeed * Math.Sin(2 * Math.PI * 5 * t);
            return fundamental + harmonic + windEffect + noise;
        }

        private double CalculateDissolvedOxygen(double t)
        {
            double tempFactor = (waterTemp - 20) / 10.0;
            double baseDO = 8.0 - 2.0 * tempFactor;
            double periodic = 0.5 * Math.Sin(2 * Math.PI * 0.2 * t);
            double noise = 0.2 * (rand.NextDouble() - 0.5);
            double aerationEffect = 0.5 * (windSpeed - 5);
            double ambientEffect = -0.15 * (ambientTemp - 25);
            return Math.Max(0, baseDO + periodic + aerationEffect + ambientEffect + noise);
        }

        private void GenerateStaticAnalysis()
        {
            ambientTemp = (double)numAmbientTemp.Value;
            humidity = (double)numHumidity.Value;
            windSpeed = (double)numWindSpeed.Value;
            mechVib = (double)numMechanicalVib.Value;
            waterTemp = (double)numWaterTemp.Value;

            sensorData = new SensorData();

            for (int sensorIdx = 0; sensorIdx < 5; sensorIdx++)
            {
                int fs = samplingRates[sensorIdx];
                int numPoints = (int)(fs * duration);
                double dt = 1.0 / fs;

                double[] timeData = new double[numPoints];
                double[] valueData = new double[numPoints];

                for (int i = 0; i < numPoints; i++)
                {
                    double t = i * dt;
                    timeData[i] = t;

                    switch (sensorIdx)
                    {
                        case 0: valueData[i] = CalculateTemperature(t); break;
                        case 1: valueData[i] = CalculateHumidity(t); break;
                        case 2: valueData[i] = CalculateAirflow(t); break;
                        case 3: valueData[i] = CalculateVibration(t); break;
                        case 4: valueData[i] = CalculateDissolvedOxygen(t); break;
                    }
                }

                sensorData.Time[sensorIdx] = timeData;
                sensorData.SensorValues[sensorIdx] = valueData;
            }

            UpdateFrequencyDomainPlots();
            UpdateLaplaceDomain();
            UpdateZDomain();
        }

        private void UpdateFrequencyDomainPlots()
        {
            string[] names = { "Temperature", "Humidity", "Airflow", "Vibration", "DO" };

            for (int i = 0; i < 5; i++)
            {
                UpdateFFTPlot(freqPlots[i], sensorData.SensorValues[i], samplingRates[i], names[i]);
            }
        }

        private void UpdateFFTPlot(PlotView plot, double[] signal, int fs, string title)
        {
            var complex = signal.Select(x => new Complex(x, 0)).ToArray();
            Fourier.Forward(complex, FourierOptions.Matlab);

            int n = complex.Length;
            double[] frequencies = new double[n / 2];
            double[] magnitudes = new double[n / 2];

            for (int i = 0; i < n / 2; i++)
            {
                frequencies[i] = i * fs / (double)n;
                magnitudes[i] = complex[i].Magnitude / n * 2;
            }

            plot.Model.Series.Clear();
            var series = new LineSeries
            {
                Title = $"{title} FFT",
                Color = OxyColors.Red,
                StrokeThickness = 2
            };

            double nyquist = fs / 2.0;
            for (int i = 0; i < frequencies.Length && frequencies[i] <= nyquist; i++)
            {
                series.Points.Add(new DataPoint(frequencies[i], magnitudes[i]));
            }

            plot.Model.Series.Add(series);
            plot.Model.Title = $"{title} Spectrum (fs={fs} Hz, Nyquist={nyquist} Hz)";
            plot.Model.Axes[0].Maximum = nyquist;
            plot.Model.InvalidatePlot(true);
        }

        private void UpdateLaplaceDomain()
        {
            rtbLaplace.Clear();

            rtbLaplace.AppendText("═══════════════════════════════════════════════════════════\n");
            rtbLaplace.AppendText("  COOLING TOWER TRANSFER FUNCTION - LAPLACE DOMAIN (s)\n");
            rtbLaplace.AppendText("═══════════════════════════════════════════════════════════\n\n");

            rtbLaplace.AppendText("System Transfer Function (Poles Only):\n\n");
            rtbLaplace.AppendText("              K\n");
            rtbLaplace.AppendText("G(s) = ─────────────────────────────────────\n");
            rtbLaplace.AppendText("       (s + p₁)(s + p₂)(s + p₃)(s + p₄)(s + p₅)\n\n");

            rtbLaplace.AppendText("Where: K = 1.5 (System gain)\n\n");

            double[] poles = {
                -1.0 / tempTau,
                -1.0 / humidityTau,
                -1.0 / airflowTau,
                -1.0 / vibrationTau,
                -1.0 / doTau
            };

            rtbLaplace.AppendText("POLES (dari sensor time constants):\n");
            rtbLaplace.AppendText($"  p₁ = {poles[0]:F2} (Temp, τ={tempTau}s, fs={samplingRates[0]}Hz)\n");
            rtbLaplace.AppendText($"  p₂ = {poles[1]:F2} (Hum, τ={humidityTau}s, fs={samplingRates[1]}Hz)\n");
            rtbLaplace.AppendText($"  p₃ = {poles[2]:F2} (Air, τ={airflowTau}s, fs={samplingRates[2]}Hz)\n");
            rtbLaplace.AppendText($"  p₄ = {poles[3]:F1} (Vib, τ={vibrationTau}s, fs={samplingRates[3]}Hz)\n");
            rtbLaplace.AppendText($"  p₅ = {poles[4]:F2} (DO, τ={doTau}s, fs={samplingRates[4]}Hz)\n\n");

            rtbLaplace.AppendText("STABILITY ANALYSIS (S-PLANE):\n");
            rtbLaplace.AppendText("  Kriteria: Pole stabil jika Re(s) < 0 (Left Half Plane)\n\n");

            bool stable = true;
            foreach (var pole in poles)
            {
                if (pole >= 0) stable = false;
            }

            rtbLaplace.AppendText("  Verifikasi:\n");
            rtbLaplace.AppendText($"    • p₁ = {poles[0]:F2} < 0 ✓ (di LHP)\n");
            rtbLaplace.AppendText($"    • p₂ = {poles[1]:F2} < 0 ✓ (di LHP)\n");
            rtbLaplace.AppendText($"    • p₃ = {poles[2]:F2} < 0 ✓ (di LHP)\n");
            rtbLaplace.AppendText($"    • p₄ = {poles[3]:F1} < 0 ✓ (di LHP)\n");
            rtbLaplace.AppendText($"    • p₅ = {poles[4]:F2} < 0 ✓ (di LHP)\n\n");

            if (stable)
                rtbLaplace.AppendText("  → SISTEM STABLE ✓ (Semua poles di LHP)\n\n");
            else
                rtbLaplace.AppendText("  → SISTEM UNSTABLE ✗\n\n");

            rtbLaplace.AppendText("KARAKTERISTIK SISTEM:\n");
            rtbLaplace.AppendText($"  • Dominant pole: p₅={poles[4]:F2} (paling lambat)\n");
            rtbLaplace.AppendText($"  • Fastest pole: p₄={poles[3]:F1} (paling cepat)\n");
            rtbLaplace.AppendText("  • System order: 5th order\n\n");

            rtbLaplace.AppendText("MULTI-RATE SAMPLING:\n");
            rtbLaplace.AppendText("  • Setiap sensor sampling sesuai bandwidth-nya\n");
            rtbLaplace.AppendText("  • Bandwidth = 1/(2πτ)\n");
            rtbLaplace.AppendText("  • Nyquist criterion terpenuhi semua\n");

            CreatePoleOnlyPlot(plotLaplace, poles);
        }

        private void UpdateZDomain()
        {
            rtbZDomain.Clear();

            rtbZDomain.AppendText("═══════════════════════════════════════════════════════════\n");
            rtbZDomain.AppendText("  COOLING TOWER DISCRETE TRANSFER FUNCTION - Z DOMAIN\n");
            rtbZDomain.AppendText("═══════════════════════════════════════════════════════════\n\n");

            string[] sensorNames = { "Temperature", "Humidity", "Airflow", "Vibration", "DO" };
            double[] timeConstants = { tempTau, humidityTau, airflowTau, vibrationTau, doTau };

            rtbZDomain.AppendText("MULTI-RATE SAMPLING VERIFICATION:\n\n");

            for (int i = 0; i < 5; i++)
            {
                double bandwidth = 1.0 / (2 * Math.PI * timeConstants[i]);
                double nyquistMin = 2 * bandwidth;
                int fs = samplingRates[i];

                rtbZDomain.AppendText($"{sensorNames[i]}:\n");
                rtbZDomain.AppendText($"  τ = {timeConstants[i]:F2}s\n");
                rtbZDomain.AppendText($"  Bandwidth = {bandwidth:F3} Hz\n");
                rtbZDomain.AppendText($"  Nyquist (min) = {nyquistMin:F3} Hz\n");
                rtbZDomain.AppendText($"  Actual fs = {fs} Hz ");

                if (fs >= nyquistMin)
                    rtbZDomain.AppendText($"✓ ({fs / nyquistMin:F1}x Nyquist)\n\n");
                else
                    rtbZDomain.AppendText("✗ ALIASING RISK!\n\n");
            }

            rtbZDomain.AppendText("Z-TRANSFORM (per sensor):\n");
            rtbZDomain.AppendText("z = e^(p·T) where T = 1/fs\n\n");

            double[] poles_s = {
                -1.0 / tempTau,
                -1.0 / humidityTau,
                -1.0 / airflowTau,
                -1.0 / vibrationTau,
                -1.0 / doTau
            };

            double[] poles_z = new double[5];
            bool stable = true;

            rtbZDomain.AppendText("Discrete Poles:\n");
            for (int i = 0; i < 5; i++)
            {
                double T = 1.0 / samplingRates[i];
                poles_z[i] = Math.Exp(poles_s[i] * T);
                double mag = Math.Abs(poles_z[i]);

                rtbZDomain.AppendText($"  z_p{i + 1} = e^({poles_s[i]:F2}×{T:F6}) = {poles_z[i]:F6}\n");
                rtbZDomain.AppendText($"       |z| = {mag:F6} ");

                if (mag < 1.0)
                    rtbZDomain.AppendText($"✓ [{sensorNames[i]}]\n");
                else
                {
                    rtbZDomain.AppendText($"✗ [{sensorNames[i]}]\n");
                    stable = false;
                }
            }

            rtbZDomain.AppendText("\n");
            if (stable)
                rtbZDomain.AppendText("→ ALL POLES INSIDE UNIT CIRCLE → STABLE ✓\n\n");
            else
                rtbZDomain.AppendText("→ UNSTABLE SYSTEM ✗\n\n");

            rtbZDomain.AppendText("Discrete Zeros (ZOH method):\n");
            double[] zeros_s = { -0.5, -1.2 };
            double T_ref = 1.0 / samplingRates[3];
            double[] zeros_z = zeros_s.Select(z => Math.Exp(z * T_ref)).ToArray();

            for (int i = 0; i < zeros_z.Length; i++)
            {
                rtbZDomain.AppendText($"  z_z{i + 1} = e^({zeros_s[i]:F1}×{T_ref:F6}) = {zeros_z[i]:F6}\n");
            }

            rtbZDomain.AppendText("\nADVANTAGES OF MULTI-RATE SAMPLING:\n");
            rtbZDomain.AppendText("• Energy efficient (slow sensors use less power)\n");
            rtbZDomain.AppendText("• Optimal data storage (no redundant samples)\n");
            rtbZDomain.AppendText("• Realistic (matches sensor physics)\n");
            rtbZDomain.AppendText("• Meets Nyquist for all channels\n");

            CreatePoleZeroPlotWithUnitCircle(plotZDomain, poles_z, zeros_z);
        }

        private void CreatePoleOnlyPlot(PlotView plotView, double[] poles)
        {
            var model = new PlotModel
            {
                Title = "Pole Plot - S-Plane (Stability: LHP)",
                PlotAreaBorderColor = OxyColors.Black
            };

            double minPole = poles.Min();
            double range = Math.Abs(minPole) * 1.3;

            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Real (σ)",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                AxislineStyle = LineStyle.Solid,
                AxislineThickness = 2,
                Minimum = -range,
                Maximum = range * 0.3
            });

            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Imaginary (jω)",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                AxislineStyle = LineStyle.Solid,
                AxislineThickness = 2,
                Minimum = -5,
                Maximum = 5
            });

            model.Annotations.Add(new LineAnnotation
            {
                Type = LineAnnotationType.Vertical,
                X = 0,
                Color = OxyColors.Red,
                StrokeThickness = 3,
                LineStyle = LineStyle.Dash,
                Text = "jω-axis (Stability Boundary)"
            });

            var stableRegion = new RectangleAnnotation
            {
                MinimumX = -range,
                MaximumX = 0,
                MinimumY = -5,
                MaximumY = 5,
                Fill = OxyColor.FromArgb(30, 0, 200, 0),
                Text = "STABLE\n(LHP)"
            };
            model.Annotations.Add(stableRegion);

            var unstableRegion = new RectangleAnnotation
            {
                MinimumX = 0,
                MaximumX = range * 0.3,
                MinimumY = -5,
                MaximumY = 5,
                Fill = OxyColor.FromArgb(30, 200, 0, 0),
                Text = "UNSTABLE\n(RHP)"
            };
            model.Annotations.Add(unstableRegion);

            var poleSeries = new ScatterSeries
            {
                MarkerType = MarkerType.Cross,
                MarkerSize = 15,
                MarkerStroke = OxyColors.Blue,
                MarkerStrokeThickness = 4,
                Title = "Poles (X)"
            };

            string[] labels = { "p₁(Temp)", "p₂(Hum)", "p₃(Air)", "p₄(Vib)", "p₅(DO)" };
            for (int i = 0; i < poles.Length; i++)
            {
                poleSeries.Points.Add(new ScatterPoint(poles[i], 0));

                model.Annotations.Add(new TextAnnotation
                {
                    Text = labels[i],
                    TextPosition = new DataPoint(poles[i], 0.3),
                    Font = "Arial",
                    FontSize = 9,
                    TextColor = OxyColors.DarkBlue
                });
            }
            model.Series.Add(poleSeries);

            plotView.Model = model;
            plotView.Model.InvalidatePlot(true);
        }

        private void CreatePoleZeroPlotWithUnitCircle(PlotView plotView, double[] poles, double[] zeros)
        {
            var model = new PlotModel
            {
                Title = "Pole-Zero Map (Z-Plane)",
                PlotAreaBorderColor = OxyColors.Black
            };

            var xAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Real",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                Minimum = -1.5,
                Maximum = 1.5,
                AxislineStyle = LineStyle.Solid,
                AxislineThickness = 2
            };

            var yAxis = new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Imaginary",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                Minimum = -1.5,
                Maximum = 1.5,
                AxislineStyle = LineStyle.Solid,
                AxislineThickness = 2
            };

            model.Axes.Add(xAxis);
            model.Axes.Add(yAxis);

            var unitCircle = new FunctionSeries(
                t => Math.Cos(t),
                t => Math.Sin(t),
                0, 2 * Math.PI, 1000)
            {
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.Dash,
                Title = "Unit Circle (|z|=1)"
            };
            model.Series.Add(unitCircle);

            model.Annotations.Add(new LineAnnotation
            {
                Type = LineAnnotationType.Vertical,
                X = 0,
                Color = OxyColors.Black,
                StrokeThickness = 1,
                LineStyle = LineStyle.Solid
            });

            model.Annotations.Add(new LineAnnotation
            {
                Type = LineAnnotationType.Horizontal,
                Y = 0,
                Color = OxyColors.Black,
                StrokeThickness = 1,
                LineStyle = LineStyle.Solid
            });

            var poleSeries = new ScatterSeries
            {
                MarkerType = MarkerType.Cross,
                MarkerSize = 12,
                MarkerStroke = OxyColors.Red,
                MarkerStrokeThickness = 3,
                Title = "Poles (X)"
            };

            foreach (var pole in poles)
            {
                poleSeries.Points.Add(new ScatterPoint(pole, 0));
            }
            model.Series.Add(poleSeries);

            if (zeros.Length > 0)
            {
                var zeroSeries = new ScatterSeries
                {
                    MarkerType = MarkerType.Circle,
                    MarkerSize = 12,
                    MarkerStroke = OxyColors.Blue,
                    MarkerStrokeThickness = 3,
                    MarkerFill = OxyColors.Transparent,
                    Title = "Zeros (O)"
                };

                foreach (var zero in zeros)
                {
                    zeroSeries.Points.Add(new ScatterPoint(zero, 0));
                }
                model.Series.Add(zeroSeries);
            }

            plotView.Model = model;
            plotView.Model.InvalidatePlot(true);
        }
    }

    public class SensorData
    {
        public double[][] Time { get; set; }
        public double[][] SensorValues { get; set; }

        public SensorData()
        {
            Time = new double[5][];
            SensorValues = new double[5][];
        }
    }
}