# 多序列比对公共网站部署说明

这个版本是公共上传网站形态：用户打开网页后上传自己的引物表 Excel 和样本 FASTA 文件，点击分析后下载完整结果压缩包。

## 启动入口

Streamlit 入口文件：

```text
app/streamlit_app.py
```

## 依赖文件

```text
requirements.txt
packages.txt
.streamlit/config.toml
```

`packages.txt` 会安装 Chromium，用于把分析生成的 HTML 图转换成 PDF。

## 推荐部署方式

### 方式 1：Streamlit Community Cloud

1. 把 `multi-sequence-aligner` 文件夹上传到 GitHub 仓库。
2. 在 Streamlit Cloud 新建应用。
3. App file 填写：

```text
app/streamlit_app.py
```

4. 部署完成后，把生成的网址发给用户。

### 方式 2：云服务器

服务器安装 Python 3.10+ 后，在本目录运行：

```bash
pip install -r requirements.txt
streamlit run app/streamlit_app.py --server.address 0.0.0.0 --server.port 8501
```

如果服务器没有 Chromium，需要安装：

```bash
sudo apt-get update
sudo apt-get install -y chromium fonts-noto-cjk
```

## 数据与隐私

上传文件和结果默认保存在服务器临时目录：

```text
系统临时目录/multi_sequence_aligner
```

如需指定保存位置，可设置环境变量：

```bash
MULTISEQ_APP_DATA=/path/to/app-data
```

公共网站建议后续增加访问控制、任务排队、定期清理临时文件和文件大小限制。
