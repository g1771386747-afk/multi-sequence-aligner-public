# 多序列比对公共网站部署说明

这个版本是公共上传网站形态：用户打开网页后上传自己的 FASTA 文件即可进行直接多序列比对；如果需要按正反引物定位理论扩增片段，也可以切换到“引物扩增片段分析”模式后上传引物表 Excel。

网站目前包含两种入口：

- 直接多序列比对：不需要引物表，适合任意同源 DNA 序列、片段序列或基因序列。支持一个多序列 FASTA，或多个 FASTA 文件。
- 引物扩增片段分析：上传引物表和样本 FASTA，生成理论扩增片段、多样本比对、可选测序验证和 PDF 结果。

直接多序列比对会自动生成：

- 输入体检表：检查空序列、重复样本名、异常碱基、长度差异和可能的反向互补序列。
- 结果概览：样本数、比对长度、变异位点、SNP、InDel、保守位点比例。
- 共识序列 FASTA：支持 IUPAC 模糊碱基。
- 样本两两相似度矩阵：CSV 表格和网页热图。
- 相对参考序列差异表：每个样本相对第一个序列的 identity、mismatch 和 gap。
- 可选 A/B 分组固定差异表。

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

网站会在用户开始新一轮分析时清理当前会话的旧结果，并自动清理超过 24 小时的会话目录。公共网站后续仍建议增加任务排队和更细的文件大小限制。若面向完全公开用户，也建议在页面上补充数据隐私提示。
