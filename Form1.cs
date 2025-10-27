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
        private Button btnSimulate, btnStopSimulation;
        private Label lblStatus;

        // Plot Views
        private readonly PlotView[] timePlots = new PlotView[5];
        private readonly PlotView[] freqPlots = new PlotView[5];
        private PlotView plotLaplace, plotZDomain;
        private RichTextBox rtbLaplace, rtbZDomain;

        // Data Storage
        private SensorData sensorData;
        private readonly int samplingRate = 1000; // Hz
        private readonly double duration = 5.0; // seconds

        // Realtime simulation
        private Timer simulationTimer;
        private double currentTime = 0;
        private readonly double timeWindow = 10.0; // seconds to display
        private bool isSimulating = false;
        private Random rand = new Random(42);

        // Parameters
        private double ambientTemp, humidity, windSpeed, mechVib, waterTemp;

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
            this.Text = "Cooling Tower Efficiency Analysis - Sensor Simulation System (REALTIME)";
            this.StartPosition = FormStartPosition.CenterScreen;
            this.ResumeLayout(false);
        }

        private void InitializeTimer()
        {
            simulationTimer = new Timer
            {
                Interval = 50 // Update every 50ms (20 Hz refresh rate)
            };
            simulationTimer.Tick += SimulationTimer_Tick;
        }

        private void InitializeUI()
        {
            // Create TabControl
            tabControl = new TabControl
            {
                Dock = DockStyle.Fill,
                Font = new Font("Segoe UI", 9F)
            };

            // Create Tabs
            tabParameters = new TabPage("Parameter Input");
            tabTimeDomain = new TabPage("Time Domain (REALTIME)");
            tabFrequencyDomain = new TabPage("Frequency Domain (FFT)");
            tabLaplace = new TabPage("S-Domain (Laplace)");
            tabZDomain = new TabPage("Z-Domain (Discrete)");

            tabControl.TabPages.AddRange(new[] { tabParameters, tabTimeDomain, tabFrequencyDomain, tabLaplace, tabZDomain });

            // Initialize each tab
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

            // Title
            Label title = new Label
            {
                Text = "Sensor Input Parameters Configuration",
                Font = new Font("Segoe UI", 14F, FontStyle.Bold),
                Location = new Point(20, yPos),
                AutoSize = true
            };
            panel.Controls.Add(title);
            yPos += 50;

            // Ambient Temperature
            AddParameterControl(panel, "Ambient Temperature (°C):", 25m, 0m, 50m, ref numAmbientTemp, ref yPos);

            // Humidity
            AddParameterControl(panel, "Relative Humidity (%):", 60m, 0m, 100m, ref numHumidity, ref yPos);

            // Wind Speed
            AddParameterControl(panel, "Wind Speed (m/s):", 5m, 0m, 30m, ref numWindSpeed, ref yPos);

            // Mechanical Vibration
            AddParameterControl(panel, "Mechanical Vibration Base (g):", 0.5m, 0m, 5m, ref numMechanicalVib, ref yPos);

            // Water Temperature
            AddParameterControl(panel, "Water Temperature (°C):", 30m, 10m, 60m, ref numWaterTemp, ref yPos);

            yPos += 20;

            // Start Simulation Button
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

            // Stop Simulation Button
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

            // Status Label
            lblStatus = new Label
            {
                Text = "Status: Stopped",
                Location = new Point(20, yPos),
                Size = new Size(400, 30),
                Font = new Font("Segoe UI", 10F, FontStyle.Bold),
                ForeColor = Color.Gray
            };
            panel.Controls.Add(lblStatus);
            yPos += 40;

            // Info Label
            Label info = new Label
            {
                Text = "Sensors: SHT85 (Temp/Humidity), Testo 440 (Airflow), PCB 352C33 (Vibration), Apera DO850 (Water Quality)\n" +
                       "Real-time simulation akan menampilkan data sensor yang terus update setiap 50ms",
                Location = new Point(20, yPos),
                Size = new Size(800, 60),
                Font = new Font("Segoe UI", 9F, FontStyle.Italic),
                ForeColor = Color.Gray
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
                "Temperature (°C) - REALTIME",
                "Humidity (%) - REALTIME",
                "Airflow (m/s) - REALTIME",
                "Vibration (g) - REALTIME",
                "Dissolved Oxygen (mg/L) - REALTIME"
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
                    Maximum = timeWindow
                });
                model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = titles[i].Split('-')[0].Trim() });

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
                "DO Spectrum"
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
                model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Frequency (Hz)" });
                model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Magnitude" });
                freqPlots[i].Model = model;

                panel.Controls.Add(freqPlots[i]);
            }

            tabFrequencyDomain.Controls.Add(panel);
        }

        private void InitializeLaplaceTab()
        {
            Panel mainPanel = new Panel { Dock = DockStyle.Fill };

            // Split: Text di kiri, grafik di kanan
            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.Fill,
                Orientation = Orientation.Vertical,
                SplitterDistance = 600
            };

            // Left: Text explanation
            rtbLaplace = new RichTextBox
            {
                Dock = DockStyle.Fill,
                Font = new Font("Consolas", 9F),
                ReadOnly = true,
                BackColor = Color.WhiteSmoke,
                Padding = new Padding(10)
            };
            split.Panel1.Controls.Add(rtbLaplace);

            // Right: Pole-Zero plot
            Panel plotPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(10) };
            Label lblPlot = new Label
            {
                Text = "Pole-Zero Map (S-Plane)",
                Dock = DockStyle.Top,
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                Height = 30,
                TextAlign = ContentAlignment.MiddleCenter
            };

            plotLaplace = new PlotView
            {
                Dock = DockStyle.Fill,
                BackColor = Color.White
            };

            plotPanel.Controls.Add(plotLaplace);
            plotPanel.Controls.Add(lblPlot);
            split.Panel2.Controls.Add(plotPanel);

            mainPanel.Controls.Add(split);
            tabLaplace.Controls.Add(mainPanel);
        }

        private void InitializeZDomainTab()
        {
            Panel mainPanel = new Panel { Dock = DockStyle.Fill };

            // Split: Text di kiri, grafik di kanan
            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.Fill,
                Orientation = Orientation.Vertical,
                SplitterDistance = 600
            };

            // Left: Text explanation
            rtbZDomain = new RichTextBox
            {
                Dock = DockStyle.Fill,
                Font = new Font("Consolas", 9F),
                ReadOnly = true,
                BackColor = Color.WhiteSmoke,
                Padding = new Padding(10)
            };
            split.Panel1.Controls.Add(rtbZDomain);

            // Right: Pole-Zero plot with unit circle
            Panel plotPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(10) };
            Label lblPlot = new Label
            {
                Text = "Pole-Zero Map (Z-Plane with Unit Circle)",
                Dock = DockStyle.Top,
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                Height = 30,
                TextAlign = ContentAlignment.MiddleCenter
            };

            plotZDomain = new PlotView
            {
                Dock = DockStyle.Fill,
                BackColor = Color.White
            };

            plotPanel.Controls.Add(plotZDomain);
            plotPanel.Controls.Add(lblPlot);
            split.Panel2.Controls.Add(plotPanel);

            mainPanel.Controls.Add(split);
            tabZDomain.Controls.Add(mainPanel);
        }

        private void BtnSimulate_Click(object sender, EventArgs e)
        {
            // Get parameters
            ambientTemp = (double)numAmbientTemp.Value;
            humidity = (double)numHumidity.Value;
            windSpeed = (double)numWindSpeed.Value;
            mechVib = (double)numMechanicalVib.Value;
            waterTemp = (double)numWaterTemp.Value;

            // Reset time
            currentTime = 0;

            // Clear all plots
            foreach (var plot in timePlots)
            {
                ((LineSeries)plot.Model.Series[0]).Points.Clear();
            }

            // Start simulation
            isSimulating = true;
            simulationTimer.Start();

            // Update UI
            btnSimulate.Enabled = false;
            btnStopSimulation.Enabled = true;
            lblStatus.Text = "Status: Running (Realtime)";
            lblStatus.ForeColor = Color.Green;

            // Also generate static analysis for FFT, Laplace, Z-domain
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

            // Increment time
            currentTime += 0.05; // 50ms

            // Calculate current sensor values
            double temp = CalculateTemperature(currentTime);
            double hum = CalculateHumidity(currentTime);
            double airflow = CalculateAirflow(currentTime);
            double vib = CalculateVibration(currentTime);
            double dox = CalculateDissolvedOxygen(currentTime);

            // Update plots
            UpdateRealtimePlot(timePlots[0], currentTime, temp);
            UpdateRealtimePlot(timePlots[1], currentTime, hum);
            UpdateRealtimePlot(timePlots[2], currentTime, airflow);
            UpdateRealtimePlot(timePlots[3], currentTime, vib);
            UpdateRealtimePlot(timePlots[4], currentTime, dox);
        }

        private void UpdateRealtimePlot(PlotView plot, double time, double value)
        {
            var series = (LineSeries)plot.Model.Series[0];
            series.Points.Add(new DataPoint(time, value));

            // Remove old points outside time window
            while (series.Points.Count > 0 && series.Points[0].X < time - timeWindow)
            {
                series.Points.RemoveAt(0);
            }

            // Update axes
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
            return ambientTemp +
                   2 * Math.Sin(2 * Math.PI * 0.5 * t) +
                   0.5 * (rand.NextDouble() - 0.5);
        }

        private double CalculateHumidity(double t)
        {
            return humidity +
                   5 * Math.Sin(2 * Math.PI * 0.3 * t) +
                   1 * (rand.NextDouble() - 0.5);
        }

        private double CalculateAirflow(double t)
        {
            return windSpeed +
                   1.5 * Math.Sin(2 * Math.PI * 1.0 * t) +
                   0.5 * Math.Sin(2 * Math.PI * 3.0 * t) +
                   0.3 * (rand.NextDouble() - 0.5);
        }

        private double CalculateVibration(double t)
        {
            return mechVib * Math.Sin(2 * Math.PI * 60 * t) +
                   0.2 * Math.Sin(2 * Math.PI * 120 * t) +
                   0.1 * (rand.NextDouble() - 0.5);
        }

        private double CalculateDissolvedOxygen(double t)
        {
            double tempFactor = (waterTemp - 20) / 10.0;
            return 8.0 - tempFactor +
                   0.5 * Math.Sin(2 * Math.PI * 0.2 * t) +
                   0.2 * (rand.NextDouble() - 0.5);
        }

        private void GenerateStaticAnalysis()
        {
            // Generate snapshot data for FFT analysis
            int numPoints = (int)(samplingRate * duration);
            double dt = 1.0 / samplingRate;

            sensorData = new SensorData
            {
                Time = new double[numPoints],
                Temperature = new double[numPoints],
                Humidity = new double[numPoints],
                Airflow = new double[numPoints],
                Vibration = new double[numPoints],
                DissolvedOxygen = new double[numPoints]
            };

            for (int i = 0; i < numPoints; i++)
            {
                double t = i * dt;
                sensorData.Time[i] = t;
                sensorData.Temperature[i] = CalculateTemperature(t);
                sensorData.Humidity[i] = CalculateHumidity(t);
                sensorData.Airflow[i] = CalculateAirflow(t);
                sensorData.Vibration[i] = CalculateVibration(t);
                sensorData.DissolvedOxygen[i] = CalculateDissolvedOxygen(t);
            }

            // Update FFT plots
            UpdateFrequencyDomainPlots();

            // Update Laplace domain
            UpdateLaplaceDomain();

            // Update Z domain
            UpdateZDomain();
        }

        private void UpdateFrequencyDomainPlots()
        {
            UpdateFFTPlot(freqPlots[0], sensorData.Temperature, "Temperature");
            UpdateFFTPlot(freqPlots[1], sensorData.Humidity, "Humidity");
            UpdateFFTPlot(freqPlots[2], sensorData.Airflow, "Airflow");
            UpdateFFTPlot(freqPlots[3], sensorData.Vibration, "Vibration");
            UpdateFFTPlot(freqPlots[4], sensorData.DissolvedOxygen, "DO");
        }

        private void UpdateFFTPlot(PlotView plot, double[] signal, string title)
        {
            // Perform FFT
            var complex = signal.Select(x => new Complex(x, 0)).ToArray();
            Fourier.Forward(complex, FourierOptions.Matlab);

            int n = complex.Length;
            double[] frequencies = new double[n / 2];
            double[] magnitudes = new double[n / 2];

            for (int i = 0; i < n / 2; i++)
            {
                frequencies[i] = i * samplingRate / (double)n;
                magnitudes[i] = complex[i].Magnitude / n * 2;
            }

            plot.Model.Series.Clear();
            var series = new LineSeries
            {
                Title = $"{title} FFT",
                Color = OxyColors.Red,
                StrokeThickness = 2
            };

            for (int i = 0; i < frequencies.Length && frequencies[i] < 100; i++)
            {
                series.Points.Add(new DataPoint(frequencies[i], magnitudes[i]));
            }

            plot.Model.Series.Add(series);
            plot.Model.InvalidatePlot(true);
        }

        private void UpdateLaplaceDomain()
        {
            // Update text
            rtbLaplace.Clear();

            rtbLaplace.AppendText("═══════════════════════════════════════════════════════════\n");
            rtbLaplace.AppendText("  COOLING TOWER TRANSFER FUNCTION - LAPLACE DOMAIN (s)\n");
            rtbLaplace.AppendText("═══════════════════════════════════════════════════════════\n\n");

            rtbLaplace.AppendText("System Transfer Function:\n\n");
            rtbLaplace.AppendText("         K(s + z₁)(s + z₂)\n");
            rtbLaplace.AppendText("G(s) = ─────────────────────────────────────\n");
            rtbLaplace.AppendText("       (s + p₁)(s + p₂)(s + p₃)(s + p₄)(s + p₅)\n\n");

            rtbLaplace.AppendText("Where:\n");
            rtbLaplace.AppendText("  K = 1.5 (System gain)\n\n");

            double[] zeros = { -0.5, -1.2 };
            rtbLaplace.AppendText("Zeros (numerator roots):\n");
            rtbLaplace.AppendText($"  z₁ = {zeros[0]:F1} (Temperature dynamics)\n");
            rtbLaplace.AppendText($"  z₂ = {zeros[1]:F1} (Humidity coupling)\n\n");

            double[] poles = { -2.0, -1.5, -5.0, -50.0, -0.8 };
            rtbLaplace.AppendText("Poles (denominator roots):\n");
            rtbLaplace.AppendText($"  p₁ = {poles[0]:F1}  (Temperature sensor)\n");
            rtbLaplace.AppendText($"  p₂ = {poles[1]:F1}  (Humidity sensor)\n");
            rtbLaplace.AppendText($"  p₃ = {poles[2]:F1}  (Airflow sensor)\n");
            rtbLaplace.AppendText($"  p₄ = {poles[3]:F1} (Vibration sensor)\n");
            rtbLaplace.AppendText($"  p₅ = {poles[4]:F1}  (Water quality sensor)\n\n");

            rtbLaplace.AppendText("System Characteristics:\n");
            rtbLaplace.AppendText("  • All poles have NEGATIVE real parts → STABLE\n");
            rtbLaplace.AppendText($"  • Dominant pole: p₅ = {poles[4]:F1} (slowest)\n");
            rtbLaplace.AppendText($"  • Fastest pole: p₄ = {poles[3]:F1} (vibration)\n");
            rtbLaplace.AppendText("  • System order: 5th order\n");

            // Update plot
            CreatePoleZeroPlot(plotLaplace, poles, zeros, "S-Plane");
        }

        private void UpdateZDomain()
        {
            double T = 1.0 / samplingRate;

            rtbZDomain.Clear();

            rtbZDomain.AppendText("═══════════════════════════════════════════════════════════\n");
            rtbZDomain.AppendText("  COOLING TOWER DISCRETE TRANSFER FUNCTION - Z DOMAIN\n");
            rtbZDomain.AppendText("═══════════════════════════════════════════════════════════\n\n");

            rtbZDomain.AppendText($"Sampling Frequency: {samplingRate} Hz\n");
            rtbZDomain.AppendText($"Sampling Period: T = {T:F6} seconds\n\n");

            rtbZDomain.AppendText("Z-Transform using Tustin (Bilinear) method:\n");
            rtbZDomain.AppendText("         2   z - 1\n");
            rtbZDomain.AppendText("    s = ─── · ─────\n");
            rtbZDomain.AppendText("         T   z + 1\n\n");

            // Calculate discrete poles
            double[] poles_s = { -2.0, -1.5, -5.0, -50.0, -0.8 };
            double[] zeros_s = { -0.5, -1.2 };

            double[] poles_z = poles_s.Select(p => Math.Exp(p * T)).ToArray();
            double[] zeros_z = zeros_s.Select(z => Math.Exp(z * T)).ToArray();

            rtbZDomain.AppendText("Discrete Poles (from s-domain):\n");
            for (int i = 0; i < poles_z.Length; i++)
            {
                rtbZDomain.AppendText($"  z_pole{i + 1} = e^({poles_s[i]:F1}·T) = {poles_z[i]:F6}\n");
            }
            rtbZDomain.AppendText("\n");

            rtbZDomain.AppendText("Discrete Zeros:\n");
            for (int i = 0; i < zeros_z.Length; i++)
            {
                rtbZDomain.AppendText($"  z_zero{i + 1} = e^({zeros_s[i]:F1}·T) = {zeros_z[i]:F6}\n");
            }
            rtbZDomain.AppendText("\n");

            rtbZDomain.AppendText("Stability Analysis:\n");
            bool stable = true;
            for (int i = 0; i < poles_z.Length; i++)
            {
                double magnitude = Math.Abs(poles_z[i]);
                rtbZDomain.AppendText($"  |z_pole{i + 1}| = {magnitude:F6} ");
                if (magnitude < 1.0)
                    rtbZDomain.AppendText("✓ (inside unit circle)\n");
                else
                {
                    rtbZDomain.AppendText("✗ (outside unit circle)\n");
                    stable = false;
                }
            }

            if (stable)
                rtbZDomain.AppendText("\n→ System Status: STABLE ✓\n");
            else
                rtbZDomain.AppendText("\n→ System Status: UNSTABLE ✗\n");

            // Update plot with unit circle
            CreatePoleZeroPlotWithUnitCircle(plotZDomain, poles_z, zeros_z);
        }

        private void CreatePoleZeroPlot(PlotView plotView, double[] poles, double[] zeros, string planeName)
        {
            var model = new PlotModel
            {
                Title = $"Pole-Zero Map ({planeName})",
                PlotAreaBorderColor = OxyColors.Black
            };

            // Axes
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Real",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray
            });

            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Imaginary",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray
            });

            // Vertical line at Re = 0 (stability boundary for s-plane)
            model.Annotations.Add(new LineAnnotation
            {
                Type = LineAnnotationType.Vertical,
                X = 0,
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.Dash,
                Text = "Stability Boundary"
            });

            // Plot poles (X markers)
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

            // Plot zeros (O markers)
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

            // Axes
            var xAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Real",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                Minimum = -1.5,
                Maximum = 1.5
            };

            var yAxis = new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Imaginary",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.LightGray,
                Minimum = -1.5,
                Maximum = 1.5
            };

            model.Axes.Add(xAxis);
            model.Axes.Add(yAxis);

            // Draw unit circle (stability boundary for z-plane)
            var unitCircle = new FunctionSeries(
                t => Math.Cos(t),
                t => Math.Sin(t),
                0, 2 * Math.PI, 1000)
            {
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.Dash,
                Title = "Unit Circle (Stability Boundary)"
            };
            model.Series.Add(unitCircle);

            // Axes lines
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

            // Plot poles (X markers)
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

            // Plot zeros (O markers)
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

    // Data structure to hold sensor readings
    public class SensorData
    {
        public double[] Time { get; set; }
        public double[] Temperature { get; set; }
        public double[] Humidity { get; set; }
        public double[] Airflow { get; set; }
        public double[] Vibration { get; set; }
        public double[] DissolvedOxygen { get; set; }
    }
}