# Sentaurus VMware SSH Workflow

本文记录在 Windows 主机上通过 VMware Workstation Pro 运行 Red Hat
Enterprise Linux 7.9 / Sentaurus 2022 虚拟机，并通过 SSH 自动上传输入、运行
仿真和取回结果的流程。

## 当前已验证配置

- Windows 仓库：`D:\code-repo\vela-tcad`
- 虚拟机操作系统：Red Hat Enterprise Linux 7.9
- Sentaurus Device：`T-2022.03-SP2`，x86_64 Linux
- Sentaurus 安装路径：`/atctools/Synopsys/tcad/T-2022.03/bin`
- ICWBEV 路径：`/atctools/Synopsys/icwbev/O-2018.06-SP2/bin`
- jEdit：`/usr/local/bin/jedit`
- SSH 别名：`sentaurus`
- SSH 用户：`tcad`
- Host-only IPv4：`192.168.119.131`
- 用户目录：`/home/tcad`
- 许可证：`27000@tcad`
- Sentaurus Workbench 数据目录：`/data`
- 自动化远端根目录：`~/sentaurus_runs/vela_oracle_2022`

已验证以下能力：

1. Windows 可以通过密钥非交互登录虚拟机。
2. 非交互 SSH 可以找到 `swb`、`sde`、`sdevice` 和 `jedit`。
3. `sdevice -h` 报告 `T-2022.03-SP2`。
4. `/data` 对 `tcad` 用户可写。
5. `tcad:27000` 许可证端口可达。
6. PN2D 0 V smoke test 成功完成，包括 SDE 和 SDevice。

## 版本与参考数据约定

仓库中的 `reference_tcad/pn2d_sentaurus2018/source` 表示历史输入 deck 和参考
数据的来源，不表示当前执行二进制仍是 Sentaurus 2018。不要因为虚拟机升级而
批量重命名或覆盖该目录。

使用 Sentaurus 2022 重新运行这些输入时：

- 输入仍可来自 `reference_tcad/pn2d_sentaurus2018/source`；
- 新结果写入 `build-release/reference_tcad/pn2d_sentaurus2022/...`；
- manifest 必须记录 `T-2022.03-SP2`；
- 只有经过审查并明确决定更新基准时，才把选定结果复制进版本控制目录。

## VMware Host-only 网络

在虚拟机关机状态下打开 `Virtual Machine Settings -> Network Adapter`：

- 勾选 `Connected`；
- 勾选 `Connect at power on`；
- 选择 `Host-only: A private network shared with the host`；
- 不选择 NAT 或 Bridged。

在 VMware Virtual Network Editor 中确认对应网络通常为 `VMnet1`，启用主机
虚拟网卡，但不配置 NAT。

RHEL 7.9 中检查接口、连接和路由：

```bash
nmcli -t -f NAME,DEVICE connection show --active
ip -4 address
ip -4 route
```

当前 Host-only 地址为 `192.168.119.131`。隔离配置的关键是路由表中没有：

```text
default via ...
```

若需要使用 NetworkManager 配置固定地址，应从 VMware 控制台操作，避免修改
连接时丢失 SSH：

```bash
nmcli connection modify "<连接名称>" \
  ipv4.method manual \
  ipv4.addresses "192.168.119.131/24" \
  ipv4.gateway "" \
  ipv4.dns "" \
  ipv4.never-default yes \
  connection.autoconnect yes
nmcli connection up "<连接名称>"
```

确认 Windows 主机可达而公网不可达：

```bash
ping -c 3 192.168.119.1
ping -c 3 8.8.8.8
```

## RHEL 7.9 SSH 服务

RHEL 7.9 使用 systemd，不再使用旧文档中的 `service` 和 `chkconfig`：

```bash
su -
systemctl enable --now sshd
systemctl status sshd
ss -tlnp | grep ':22'
```

如果 firewalld 正在运行，允许 SSH 服务：

```bash
firewall-cmd --permanent --add-service=ssh
firewall-cmd --reload
```

## Windows SSH 配置

配置文件：

```text
C:\Users\qzw\.ssh\config
```

当前别名配置应包含：

```sshconfig
Host sentaurus
    HostName 192.168.119.131
    User tcad
    IdentityFile C:/Users/qzw/.ssh/id_rsa_sentaurus
    IdentitiesOnly yes
    ServerAliveInterval 30
    ServerAliveCountMax 3
```

`IdentityFile` 应与实际使用的私钥一致。RHEL 7.9 当前连接不需要默认放开
`ssh-rsa` 或 `ssh-dss` 弱算法；只有出现明确的算法协商错误时，才为该 Host
添加最小范围的兼容设置。

如果重装虚拟机后复用了同一 IP，应清理旧主机指纹并重新核对新指纹：

```powershell
ssh-keygen -R 192.168.119.131
ssh sentaurus
```

SSH config 应保存为 ASCII 或 UTF-8 无 BOM。若 OpenSSH 报
`no argument after keyword "\377\376"`，说明文件可能被保存成 UTF-16。

## SSH 密钥权限

虚拟机中 `tcad` 用户的权限应为：

```bash
mkdir -p /home/tcad/.ssh
touch /home/tcad/.ssh/authorized_keys
chown -R tcad:tcad /home/tcad/.ssh
chmod go-w /home/tcad
chmod 700 /home/tcad/.ssh
chmod 600 /home/tcad/.ssh/authorized_keys
```

Windows 连通性验证：

```powershell
ssh sentaurus 'hostname; pwd; whoami'
```

期望用户和目录分别为 `tcad`、`/home/tcad`。

## Sentaurus 2022 环境

为了同时支持登录 shell 和 `ssh sentaurus 'command'` 非交互 shell，建议把
Sentaurus 环境集中在：

```text
/etc/profile.d/sentaurus2022.sh
```

内容为：

```bash
export SNPSLMD_LICENSE_FILE="27000@tcad"
export ICWBEV_USER="SENTAURUS"
export STDB="/data"
export JAVA_HOME="/usr/lib/jvm/java-11-openjdk-11.0.8.10-1.el7.x86_64"
export PATH="/atctools/Synopsys/tcad/T-2022.03/bin:/atctools/Synopsys/icwbev/O-2018.06-SP2/bin:$JAVA_HOME/bin:$PATH"
```

设置系统文件权限：

```bash
chown root:root /etc/profile.d/sentaurus2022.sh
chmod 644 /etc/profile.d/sentaurus2022.sh
```

如果这些变量原先直接写在 `/etc/profile`，应移入上述独立文件，避免为了
非交互运行而重复加载整个 `/etc/profile`。

在 `/home/tcad/.bashrc` 最前面、任何非交互 `return` 之前加入：

```bash
[ -r /etc/profile.d/sentaurus2022.sh ] && . /etc/profile.d/sentaurus2022.sh
```

这一步是自动 runner 的必要条件。仅在交互登录后能够运行 `sdevice` 不够；
下列非交互检查也必须成功：

```powershell
ssh sentaurus 'command -v swb; command -v sde; command -v sdevice; command -v jedit'
ssh sentaurus 'echo "$SNPSLMD_LICENSE_FILE"; echo "$ICWBEV_USER"; echo "$STDB"; echo "$JAVA_HOME"'
ssh sentaurus 'sdevice -h 2>&1 | sed -n "1,5p"'
```

当前已验证的 `sdevice` banner 包含：

```text
Sentaurus Device
Version T-2022.03-SP2
(0.7745337, x86_64, Linux)
```

`STROOT` 不是当前安装配置的一部分；只要上述 PATH 和命令检查成功，
`STROOT` 为空不表示错误。

### Java 路径

系统升级 Java 后，版本化的 `JAVA_HOME` 可能改变。使用以下命令核对：

```bash
readlink -f /usr/bin/java
test -d "$JAVA_HOME" && echo "JAVA_HOME OK"
```

### STDB 权限

```powershell
ssh sentaurus 'test -w "$STDB" && echo "STDB writable" || echo "STDB not writable"'
```

当前已验证输出：

```text
STDB writable
```

### 许可证连通性

系统未提供 `lmutil`，使用 Bash TCP 检查许可证端口：

```powershell
ssh sentaurus "timeout 3 bash -c '</dev/tcp/tcad/27000' && echo 'license port reachable' || echo 'license port unreachable'"
```

当前已验证输出：

```text
license port reachable
```

端口可达只验证网络；真实 SDevice smoke test 成功才证明所需许可证 feature
能够 checkout。

## PN2D 0 V smoke test

先在 Windows PowerShell 中进入仓库：

```powershell
Set-Location D:\code-repo\vela-tcad
```

先生成计划，不连接虚拟机：

```powershell
& D:\msys64\ucrt64\bin\python.exe scripts\run_sentaurus_vm_reference.py pn2d `
  --ssh-target sentaurus `
  --sentaurus-version T-2022.03-SP2 `
  --stages 0v `
  --remote-root ~/sentaurus_runs/vela_oracle_2022 `
  --local-output-dir build-release\reference_tcad\pn2d_sentaurus2022\sentaurus_vm_runs `
  --run-id sentaurus2022_license_smoke `
  --dry-run
```

审查 manifest 后，去掉 `--dry-run` 执行真实 smoke test：

```powershell
& D:\msys64\ucrt64\bin\python.exe scripts\run_sentaurus_vm_reference.py pn2d `
  --ssh-target sentaurus `
  --sentaurus-version T-2022.03-SP2 `
  --stages 0v `
  --remote-root ~/sentaurus_runs/vela_oracle_2022 `
  --local-output-dir build-release\reference_tcad\pn2d_sentaurus2022\sentaurus_vm_runs `
  --run-id sentaurus2022_license_smoke
```

该 smoke test 已在当前 RHEL 7.9 / T-2022.03-SP2 虚拟机上成功完成。

结果保存在：

```text
build-release/reference_tcad/pn2d_sentaurus2022/sentaurus_vm_runs/
  sentaurus2022_license_smoke/source/
```

远端工作目录为：

```text
~/sentaurus_runs/vela_oracle_2022/sentaurus2022_license_smoke/source
```

## 手动上传和运行

历史输入目录：

```text
D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source
```

使用带时间戳的隔离目录：

```powershell
$local = "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$remote = "~/sentaurus_runs/vela_oracle_2022/manual_$stamp/source"

ssh sentaurus "mkdir -p $remote"
scp -r "$local\*" "sentaurus:$remote/"

ssh sentaurus "cd $remote && sde -e -l pn2d_sde.cmd > run_sde.out 2>&1"
ssh sentaurus "cd $remote && sdevice pn2d_0v_sdevice.cmd > run_pn2d_0v.out 2>&1"
ssh sentaurus "cd $remote && tail -n 80 run_pn2d_0v.out"
```

不要把生成文件直接取回到
`reference_tcad/pn2d_sentaurus2018/source`。应先复制到
`build-release/reference_tcad/pn2d_sentaurus2022/...` 中审查。

## 常见问题

### 手工登录能运行，远程命令找不到 sdevice

症状：

```powershell
ssh sentaurus 'command -v sdevice'
```

没有输出，但手工 `ssh sentaurus` 后能够运行。

原因是 Sentaurus 环境只由交互或登录 shell 加载。确认
`/home/tcad/.bashrc` 在任何早期 `return` 之前加载
`/etc/profile.d/sentaurus2022.sh`。

### sdevice 可以显示帮助，但仿真许可证失败

依次检查：

```powershell
ssh sentaurus 'echo "$SNPSLMD_LICENSE_FILE"'
ssh sentaurus "timeout 3 bash -c '</dev/tcp/tcad/27000' && echo reachable || echo unreachable"
```

随后检查真实运行日志中的 `license`、`checkout`、`error` 和 `fatal`：

```powershell
ssh sentaurus 'grep -Ei "license|checkout|error|fatal" ~/sentaurus_runs/vela_oracle_2022/<run-id>/source/run_pn2d_0v.out'
```

### 虚拟机突然连不上

从 VMware 控制台检查：

```bash
ip -4 address
ip -4 route
systemctl status sshd
ss -tlnp | grep ':22'
```

若 IP 发生变化，更新 `C:\Users\qzw\.ssh\config` 中的 `HostName`，以及文档和
任何显式传入 `--host-name` 的专用脚本参数。

## Vela Runner Contract

自动 runner 当前约定：

- SSH target：`sentaurus`
- 当前执行版本：通过 `--sentaurus-version T-2022.03-SP2` 写入 manifest
- Remote root：`~/sentaurus_runs/vela_oracle_2022`
- Local PN2D input：`reference_tcad/pn2d_sentaurus2018/source`
- Sentaurus 2022 staged output：
  `build-release/reference_tcad/pn2d_sentaurus2022/sentaurus_vm_runs`
- Required common files：`pn2d_sde.cmd`、`models.par`
- 0 V command：
  `sdevice pn2d_0v_sdevice.cmd > run_pn2d_0v.out 2>&1`

runner 绝不能覆盖 `reference_tcad/pn2d_sentaurus2018/source`。所有新产物先进入
`build-release` 下的版本化 staging 目录，经审查后再决定是否更新参考基准。
