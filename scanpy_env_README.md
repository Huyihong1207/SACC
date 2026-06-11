环境位置：`/Users/niubi/Desktop/SACC/.venv-scanpy`

激活：

```bash
source /Users/niubi/Desktop/SACC/.venv-scanpy/bin/activate
```

验证：

```bash
python -c "import scanpy as sc, harmonypy; print(sc.__version__, harmonypy.__version__)"
```

说明：

- `git` 已存在于系统中：`/usr/bin/git`
- 当前环境已按单细胞分析常用依赖配置，包含 `scanpy`、`anndata`、`harmonypy`、`igraph`、`leidenalg`、`jupyterlab`
- 未安装 `tables`：该包需要本机 HDF5 开发头文件，当前机器缺失；通常不影响常规 `scanpy`/`h5ad` 工作流
