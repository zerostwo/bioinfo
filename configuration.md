---
title: "生信分析服务器配置"
author: "duansq"
output:
  html_document:
    toc: yes
---

拿到一台新的服务器生信分析环境的配置

# 1 终端美化（Oh My Zsh）

>参考[Ubuntu 下 Oh My Zsh 的最佳实践「安装及配置」](https://segmentfault.com/a/1190000015283092)

## 1.1 安装 Zsh

```bash
sudo apt install zsh
chsh -s /bin/zsh
```

## 1.2 安装 Oh My Zsh

```bash
sudo apt install git
wget https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | sh
```

## 1.3 插件配置

- zsh-autosuggestions(命令行命令键入时的历史命令建议插件)
  ```bash
  git clone https://github.com/zsh-users/zsh-autosuggestions${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-autosuggestions
  ```

- zsh-syntax-highlighting(命令行语法高亮插件)
  ```bash
  git clone https://github.com/zsh-users/zsh-syntax-highlighting.git${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
  ```

## 1.4 `.zshrc`配置

```vim
# 设置命令行的主题
ZSH_THEME="gianu"

# 启动错误命令自动更正
ENABLE_CORRECTION="true"

# 在命令执行的过程中，使用小红点进行提示
COMPLETION_WAITING_DOTS="true"

# 启用已安装的插件
plugins=(
  git zsh-autosuggestions zsh-syntax-highlighting
)
```
修改完`.zshrc`后记得：
```bash
source .zshrc
```
# 2 R及Rstudio的安装

## 2.1 R安装

默认安装的不是最新版本的，所以为来了获取最新版本的按照如下操作：

- 添加密钥
  ```bash
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys   E298A3A825C0D65DFD57CBB651716619E084DAB9
  ```
  
- 添加源，下面这个源是3.6的，如果需要获取其它版本请访问[UBUNTU PACKAGES FOR R](https://cran.r-project.org/bin/linux/ubuntu/README.html)
  ```bash
  sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran35/"
  ```
  
- 安装
  ```bash
  sudo apt update
  sudo apt intall r-base
  ```

## 2.2 安装Rstudio

>参考[Download RStudio Server for Debian & Ubuntu](https://rstudio.com/products/rstudio/download-server/debian-ubuntu/)

**注意：**注意不同Ubuntu系统版本下载的Rstudio Server不同。

现在我常用Ubuntu 18.04，因此下面参考流程以此版本为例：

```bash
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.5042-amd64.deb
sudo gdebi rstudio-server-1.2.5042-amd64.deb
```

Rstudio Server默认端口`8787`




