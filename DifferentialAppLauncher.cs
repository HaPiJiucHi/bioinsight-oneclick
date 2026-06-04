using System;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Net;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace DifferentialAnalysisLauncher
{
    internal static class Program
    {
        [STAThread]
        private static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new LauncherForm());
        }
    }

    public sealed class LauncherForm : Form
    {
        private readonly Label statusLabel;
        private readonly TextBox logBox;
        private readonly Button startButton;
        private readonly Button openButton;
        private readonly Button stopButton;
        private readonly Button depsButton;
        private Process rProcess;
        private int port = 3838;
        private string appDirectory;
        private string rscriptPath;
        private const string RInstallerUrl = "https://cran.r-project.org/bin/windows/base/R-latest-win.exe";

        public LauncherForm()
        {
            Text = "差异分析软件";
            StartPosition = FormStartPosition.CenterScreen;
            Size = new Size(720, 460);
            MinimumSize = new Size(640, 390);
            Font = new Font("Microsoft YaHei UI", 9F);

            var title = new Label
            {
                AutoSize = false,
                Text = "差异分析一键软件",
                Font = new Font("Microsoft YaHei UI", 18F, FontStyle.Bold),
                Location = new Point(22, 18),
                Size = new Size(520, 42)
            };

            statusLabel = new Label
            {
                AutoSize = false,
                Text = "状态：未启动",
                Location = new Point(24, 66),
                Size = new Size(650, 28)
            };

            startButton = new Button
            {
                Text = "启动软件",
                Location = new Point(24, 108),
                Size = new Size(124, 36)
            };
            startButton.Click += async (sender, args) => await StartAppAsync();

            openButton = new Button
            {
                Text = "打开界面",
                Location = new Point(162, 108),
                Size = new Size(124, 36),
                Enabled = false
            };
            openButton.Click += (sender, args) => OpenBrowser();

            stopButton = new Button
            {
                Text = "停止服务",
                Location = new Point(300, 108),
                Size = new Size(124, 36),
                Enabled = false
            };
            stopButton.Click += (sender, args) => StopServer();

            depsButton = new Button
            {
                Text = "检查依赖",
                Location = new Point(438, 108),
                Size = new Size(124, 36)
            };
            depsButton.Click += async (sender, args) => await InstallDependenciesAsync();

            logBox = new TextBox
            {
                Location = new Point(24, 164),
                Size = new Size(656, 230),
                Multiline = true,
                ScrollBars = ScrollBars.Vertical,
                ReadOnly = true,
                BackColor = Color.White
            };

            Controls.Add(title);
            Controls.Add(statusLabel);
            Controls.Add(startButton);
            Controls.Add(openButton);
            Controls.Add(stopButton);
            Controls.Add(depsButton);
            Controls.Add(logBox);

            FormClosing += (sender, args) => StopServer();
            Shown += async (sender, args) => await StartAppAsync();

            appDirectory = AppDomain.CurrentDomain.BaseDirectory.TrimEnd(Path.DirectorySeparatorChar);
            rscriptPath = FindRscript();
            AppendLog("程序目录：" + appDirectory);
            AppendLog(rscriptPath == null ? "未找到 Rscript.exe。" : "Rscript：" + rscriptPath);
        }

        private async Task StartAppAsync()
        {
            startButton.Enabled = false;
            depsButton.Enabled = false;

            try
            {
                if (!File.Exists(Path.Combine(appDirectory, "app.R")))
                {
                    throw new FileNotFoundException("找不到 app.R，请把 exe 放在 pipeline 软件目录中运行。");
                }

                rscriptPath = FindRscript();
                if (rscriptPath == null)
                {
                    throw new FileNotFoundException("找不到 Rscript.exe。请点击“检查依赖”，软件会下载 R 并安装到当前软件目录。");
                }

                port = await FindFreePortAsync(3838, 3858);
                if (await IsServerReadyAsync(port))
                {
                    SetStatus("状态：服务已在运行，正在打开界面");
                    openButton.Enabled = true;
                    startButton.Enabled = true;
                    depsButton.Enabled = true;
                    stopButton.Enabled = false;
                    OpenBrowser();
                    return;
                }

                string appPathForR = appDirectory.Replace("\\", "/");
                string rExpression = string.Format(
                    "shiny::runApp('{0}', launch.browser = FALSE, host = '127.0.0.1', port = {1})",
                    appPathForR.Replace("'", "\\'"),
                    port
                );

                var startInfo = new ProcessStartInfo
                {
                    FileName = rscriptPath,
                    Arguments = "-e \"" + rExpression.Replace("\"", "\\\"") + "\"",
                    WorkingDirectory = appDirectory,
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true
                };

                rProcess = new Process { StartInfo = startInfo, EnableRaisingEvents = true };
                rProcess.OutputDataReceived += (sender, args) => AppendLogFromThread(args.Data);
                rProcess.ErrorDataReceived += (sender, args) => AppendLogFromThread(args.Data);
                rProcess.Exited += (sender, args) => BeginInvoke(new Action(() =>
                {
                    SetStatus("状态：服务已停止");
                    openButton.Enabled = false;
                    stopButton.Enabled = false;
                    startButton.Enabled = true;
                    depsButton.Enabled = true;
                }));

                AppendLog("正在启动 Shiny 服务...");
                rProcess.Start();
                rProcess.BeginOutputReadLine();
                rProcess.BeginErrorReadLine();

                bool ready = await WaitForServerAsync(port, TimeSpan.FromSeconds(45));
                if (!ready)
                {
                    throw new TimeoutException("服务启动超时。请点击“检查依赖”，或查看日志中的 R 报错。");
                }

                SetStatus("状态：运行中  http://127.0.0.1:" + port);
                openButton.Enabled = true;
                stopButton.Enabled = true;
                OpenBrowser();
            }
            catch (Exception ex)
            {
                AppendLog("启动失败：" + ex.Message);
                SetStatus("状态：启动失败");
                startButton.Enabled = true;
                depsButton.Enabled = true;
            }
        }

        private async Task InstallDependenciesAsync()
        {
            depsButton.Enabled = false;
            startButton.Enabled = false;
            try
            {
                rscriptPath = FindRscript();
                if (rscriptPath == null)
                {
                    await InstallLocalRAsync();
                    rscriptPath = FindRscript();
                    if (rscriptPath == null)
                    {
                        throw new FileNotFoundException("R 安装后仍找不到 Rscript.exe。");
                    }
                }
                string installer = Path.Combine(appDirectory, "install_dependencies.R");
                if (!File.Exists(installer))
                {
                    throw new FileNotFoundException("找不到 install_dependencies.R。");
                }

                SetStatus("状态：正在检查依赖");
                AppendLog("开始检查/安装 R 依赖...");
                await RunProcessAsync(rscriptPath, "\"" + installer + "\"", appDirectory);
                SetStatus("状态：依赖检查完成");
                AppendLog("依赖检查完成。");
            }
            catch (Exception ex)
            {
                AppendLog("依赖检查失败：" + ex.Message);
                SetStatus("状态：依赖检查失败");
            }
            finally
            {
                depsButton.Enabled = true;
                startButton.Enabled = true;
            }
        }

        private async Task InstallLocalRAsync()
        {
            string localRDir = Path.Combine(appDirectory, "R");
            string downloadsDir = Path.Combine(appDirectory, "downloads");
            Directory.CreateDirectory(localRDir);
            Directory.CreateDirectory(downloadsDir);

            string installerPath = Path.Combine(downloadsDir, "R-latest-win.exe");
            SetStatus("状态：正在下载 R");
            AppendLog("未检测到 R，开始下载：" + RInstallerUrl);

            await Task.Run(() =>
            {
                using (var client = new WebClient())
                {
                    client.DownloadFile(RInstallerUrl, installerPath);
                }
            });

            SetStatus("状态：正在安装 R 到软件目录");
            AppendLog("正在安装 R 到：" + localRDir);
            string args = "/VERYSILENT /SUPPRESSMSGBOXES /NORESTART /DIR=\"" + localRDir + "\"";
            await RunProcessAsync(installerPath, args, appDirectory);
            AppendLog("R 安装完成。");
        }

        private Task RunProcessAsync(string fileName, string arguments, string workingDirectory)
        {
            var tcs = new TaskCompletionSource<bool>();
            var process = new Process
            {
                StartInfo = new ProcessStartInfo
                {
                    FileName = fileName,
                    Arguments = arguments,
                    WorkingDirectory = workingDirectory,
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true
                },
                EnableRaisingEvents = true
            };

            process.OutputDataReceived += (sender, args) => AppendLogFromThread(args.Data);
            process.ErrorDataReceived += (sender, args) => AppendLogFromThread(args.Data);
            process.Exited += (sender, args) =>
            {
                int exitCode = process.ExitCode;
                process.Dispose();
                if (exitCode == 0)
                {
                    tcs.TrySetResult(true);
                }
                else
                {
                    tcs.TrySetException(new InvalidOperationException("R 进程退出码：" + exitCode));
                }
            };

            process.Start();
            process.BeginOutputReadLine();
            process.BeginErrorReadLine();
            return tcs.Task;
        }

        private async Task<int> FindFreePortAsync(int start, int end)
        {
            for (int candidate = start; candidate <= end; candidate++)
            {
                if (await IsServerReadyAsync(candidate))
                {
                    return candidate;
                }
                try
                {
                    var listener = new System.Net.Sockets.TcpListener(IPAddress.Loopback, candidate);
                    listener.Start();
                    listener.Stop();
                    return candidate;
                }
                catch
                {
                }
            }
            return start;
        }

        private async Task<bool> WaitForServerAsync(int serverPort, TimeSpan timeout)
        {
            DateTime deadline = DateTime.Now.Add(timeout);
            while (DateTime.Now < deadline)
            {
                if (await IsServerReadyAsync(serverPort))
                {
                    return true;
                }
                await Task.Delay(750);
            }
            return false;
        }

        private Task<bool> IsServerReadyAsync(int serverPort)
        {
            return Task.Run(() =>
            {
                try
                {
                    var request = (HttpWebRequest)WebRequest.Create("http://127.0.0.1:" + serverPort);
                    request.Timeout = 800;
                    request.ReadWriteTimeout = 800;
                    using (var response = (HttpWebResponse)request.GetResponse())
                    {
                        return (int)response.StatusCode >= 200 && (int)response.StatusCode < 500;
                    }
                }
                catch
                {
                    return false;
                }
            });
        }

        private void OpenBrowser()
        {
            string url = "http://127.0.0.1:" + port;
            try
            {
                Process.Start(new ProcessStartInfo
                {
                    FileName = url,
                    UseShellExecute = true
                });
            }
            catch (Exception ex)
            {
                AppendLog("无法打开浏览器：" + ex.Message);
            }
        }

        private void StopServer()
        {
            try
            {
                if (rProcess != null && !rProcess.HasExited)
                {
                    AppendLog("正在停止服务...");
                    rProcess.Kill();
                    rProcess.WaitForExit(3000);
                }
            }
            catch (Exception ex)
            {
                AppendLog("停止服务失败：" + ex.Message);
            }
            finally
            {
                rProcess = null;
                openButton.Enabled = false;
                stopButton.Enabled = false;
                startButton.Enabled = true;
                depsButton.Enabled = true;
            }
        }

        private string FindRscript()
        {
            if (!string.IsNullOrEmpty(appDirectory))
            {
                string localRRoot = Path.Combine(appDirectory, "R");
                string[] localCandidates = new[]
                {
                    Path.Combine(localRRoot, "bin", "Rscript.exe"),
                    Path.Combine(localRRoot, "bin", "x64", "Rscript.exe")
                };

                foreach (string candidate in localCandidates)
                {
                    if (File.Exists(candidate))
                    {
                        return candidate;
                    }
                }

                if (Directory.Exists(localRRoot))
                {
                    foreach (string file in Directory.GetFiles(localRRoot, "Rscript.exe", SearchOption.AllDirectories))
                    {
                        return file;
                    }
                }
            }

            string[] candidates = new[]
            {
                Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles), "R", "R-4.5.3", "bin", "Rscript.exe"),
                Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles), "R", "R-4.5.3", "bin", "x64", "Rscript.exe")
            };

            foreach (string candidate in candidates)
            {
                if (File.Exists(candidate))
                {
                    return candidate;
                }
            }

            string programFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            string rRoot = Path.Combine(programFiles, "R");
            if (Directory.Exists(rRoot))
            {
                foreach (string file in Directory.GetFiles(rRoot, "Rscript.exe", SearchOption.AllDirectories))
                {
                    return file;
                }
            }

            string path = Environment.GetEnvironmentVariable("PATH") ?? "";
            foreach (string dir in path.Split(Path.PathSeparator))
            {
                try
                {
                    string candidate = Path.Combine(dir.Trim(), "Rscript.exe");
                    if (File.Exists(candidate))
                    {
                        return candidate;
                    }
                }
                catch
                {
                }
            }

            return null;
        }

        private void SetStatus(string text)
        {
            statusLabel.Text = text;
        }

        private void AppendLogFromThread(string text)
        {
            if (string.IsNullOrWhiteSpace(text))
            {
                return;
            }
            BeginInvoke(new Action(() => AppendLog(text)));
        }

        private void AppendLog(string text)
        {
            if (string.IsNullOrWhiteSpace(text))
            {
                return;
            }
            logBox.AppendText("[" + DateTime.Now.ToString("HH:mm:ss") + "] " + text + Environment.NewLine);
        }
    }
}
